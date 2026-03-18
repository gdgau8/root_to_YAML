#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <limits>
#include <map>

using namespace std;

// ============================================================
// Helpers
// ============================================================

string Trim(const string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

bool StartsWith(const string &s, const string &prefix) {
    return s.rfind(prefix, 0) == 0;
}

vector<string> SplitWS(const string &line) {
    vector<string> tokens;
    istringstream iss(line);
    string tok;
    while (iss >> tok) tokens.push_back(tok);
    return tokens;
}

bool FileExists(const string &fname) {
    ifstream f(fname.c_str());
    return f.good();
}

bool IsNanString(string s) {
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    return (s == "nan" || s == ".nan" || s == "+nan" || s == "-nan");
}

bool IsNumericToken(const string &s) {
    if (s.empty()) return false;
    if (IsNanString(s)) return true;
    char *endptr = nullptr;
    strtod(s.c_str(), &endptr);
    return (*endptr == '\0');
}

long double StringToLD(const string &s) {
    if (IsNanString(s)) return numeric_limits<long double>::quiet_NaN();
    return stold(s);
}

string NormalizeNumberString(const string &s) {
    if (IsNanString(s)) return ".nan";

    string out = s;
    if (out.find('e') == string::npos && out.find('E') == string::npos &&
        out.find('.') != string::npos) {
        while (!out.empty() && out.back() == '0') out.pop_back();
        if (!out.empty() && out.back() == '.') out.pop_back();
    }
    if (out.empty()) out = "0";
    return out;
}

string FormatFixed(long double x, int decimals) {
    ostringstream oss;
    oss << fixed << setprecision(decimals) << (double)x;
    return NormalizeNumberString(oss.str());
}

long double RoundHalfAwayFromZero(long double x) {
    if (x >= 0.0L) return floorl(x + 0.5L);
    return ceill(x - 0.5L);
}

string RoundToPlaces(const string &s, int decimals) {
    if (IsNanString(s)) return ".nan";

    long double x = StringToLD(s);
    if (!isfinite((double)x)) return NormalizeNumberString(s);

    long double factor = powl(10.0L, decimals);
    long double y = x * factor;
    long double yr = RoundHalfAwayFromZero(y);
    long double xr = yr / factor;

    return FormatFixed(xr, max(0, decimals));
}

string YAMLQuote(const string &s) {
    string out = "\"";
    for (char c : s) {
        if (c == '\\' || c == '"') out += '\\';
        out += c;
    }
    out += "\"";
    return out;
}

// ============================================================
// Data structures
// ============================================================

struct TH1VBlock {
    string raw_header_line;
    string hist_name;
    int nbins_declared = -1;
    vector<vector<string>> data; // xlow xhigh value error
};

struct FormattedRow {
    string xlow;
    string xhigh;
    string value;
    string error;
};

// ============================================================
// Parse exact input format
// ============================================================

vector<TH1VBlock> ReadTH1VBlocks(const string &file_name) {
    ifstream fin(file_name.c_str());
    vector<TH1VBlock> blocks;
    string line;

    TH1VBlock current;
    bool inBlock = false;

    auto finalizeBlock = [&]() {
        if (inBlock) {
            if (!current.hist_name.empty() && !current.data.empty()) {
                blocks.push_back(current);
            }
            current = TH1VBlock();
            inBlock = false;
        }
    };

    while (getline(fin, line)) {
        string t = Trim(line);
        if (t.empty()) continue;

        if (StartsWith(t, "#")) {
            finalizeBlock();

            vector<string> toks = SplitWS(t);
            if (toks.size() >= 4 && toks[0] == "#" && toks[1] == "TH1V") {
                current.raw_header_line = t;
                current.hist_name = toks[2];
                try {
                    current.nbins_declared = stoi(toks[3]);
                } catch (...) {
                    current.nbins_declared = -1;
                }
                inBlock = true;
            }
        } else {
            if (!inBlock) continue;

            vector<string> toks = SplitWS(t);
            if (toks.size() != 4) continue;

            bool ok = true;
            for (const auto &tok : toks) {
                if (!IsNumericToken(tok)) {
                    ok = false;
                    break;
                }
            }
            if (ok) current.data.push_back(toks);
        }
    }

    finalizeBlock();
    return blocks;
}

// ============================================================
// PDG-style formatting
// ============================================================

struct ErrResult {
    string realerr;
    int sigfigs;
    int decimals;
};

ErrResult ERR_Format(const string &err_str) {
    ErrResult res;

    if (IsNanString(err_str)) {
        res.realerr = ".nan";
        res.sigfigs = 1;
        res.decimals = 0;
        return res;
    }

    long double err = StringToLD(err_str);

    if (!isfinite((double)err) || err <= 0.0L) {
        res.realerr = NormalizeNumberString(err_str);
        res.sigfigs = 1;
        res.decimals = 0;
        return res;
    }

    int exponent = (int)floor(log10((double)fabsl(err)));
    long double first3_ld = err * powl(10.0L, -exponent + 2);
    int first3 = (int)floor(first3_ld + 1e-12L);

    int sigfigs;
    if (first3 < 355) sigfigs = 2;
    else if (first3 < 950) sigfigs = 1;
    else sigfigs = 1;

    int decimals = sigfigs - 1 - exponent;

    res.realerr = RoundToPlaces(err_str, decimals);
    res.sigfigs = sigfigs;
    res.decimals = decimals;
    return res;
}

vector<FormattedRow> FormatTH1VBlock(const TH1VBlock &block) {
    vector<FormattedRow> out;

    for (const auto &row : block.data) {
        string xlow = NormalizeNumberString(row[0]);
        string xhigh = NormalizeNumberString(row[1]);
        string value = row[2];
        string err = row[3];

        ErrResult er = ERR_Format(err);

        string rounded_value;
        if (IsNanString(err)) {
            rounded_value = NormalizeNumberString(value);
        } else if (!isfinite((double)StringToLD(err)) || StringToLD(err) <= 0.0L) {
            rounded_value = NormalizeNumberString(value);
        } else {
            rounded_value = RoundToPlaces(value, er.decimals);
        }

        FormattedRow fr;
        fr.xlow = xlow;
        fr.xhigh = xhigh;
        fr.value = rounded_value;
        fr.error = er.realerr;
        out.push_back(fr);
    }

    return out;
}

// ============================================================
// Group by identical binning
// ============================================================

string MakeBinningKey(const TH1VBlock &block) {
    ostringstream oss;
    for (const auto &row : block.data) {
        oss << NormalizeNumberString(row[0]) << ":" << NormalizeNumberString(row[1]) << "|";
    }
    return oss.str();
}

// ============================================================
// Write one combined HEPData document per group
// ============================================================

bool WriteCombinedHEPDataDocuments(
    const string &outFile,
    const vector<vector<TH1VBlock>> &groupedBlocks,
    const vector<vector<vector<FormattedRow>>> &groupedFormatted,
    const string &xName,
    const string &xUnits,
    const string &yUnits,
    const string &errorLabel,
    const string &qualifierName,
    const string &qualifierValue,
    const string &descriptionPrefix,
    const string &tableNamePrefix
) {
    ofstream fout(outFile.c_str());
    if (!fout.is_open()) return false;

    for (size_t g = 0; g < groupedBlocks.size(); g++) {
        fout << "---\n";

        // dependent_variables first
        fout << "dependent_variables:\n";
        for (size_t i = 0; i < groupedBlocks[g].size(); i++) {
            fout << "- header: {name: " << YAMLQuote(groupedBlocks[g][i].hist_name);
            if (!yUnits.empty()) fout << ", units: " << YAMLQuote(yUnits);
            fout << "}\n";

            if (!qualifierName.empty() && !qualifierValue.empty()) {
                fout << "  qualifiers:\n";
                fout << "  - {name: " << YAMLQuote(qualifierName)
                     << ", value: " << YAMLQuote(qualifierValue) << "}\n";
            }

            fout << "  values:\n";
            for (const auto &row : groupedFormatted[g][i]) {
                fout << "  - errors:\n";
                fout << "    - {symerror: " << row.error << ", label: "
                     << YAMLQuote(errorLabel) << "}\n";
                fout << "    value: " << row.value << "\n";
            }
        }

        // description
        string desc;
        if (groupedBlocks[g].size() == 1) {
            desc = descriptionPrefix + " " + groupedBlocks[g][0].hist_name + ".";
        } else {
            desc = descriptionPrefix + " ";
            for (size_t i = 0; i < groupedBlocks[g].size(); i++) {
                desc += groupedBlocks[g][i].hist_name;
                if (i + 2 < groupedBlocks[g].size()) desc += ", ";
                else if (i + 1 < groupedBlocks[g].size()) desc += ", and ";
            }
            desc += ".";
        }
        fout << "description: " << YAMLQuote(desc) << "\n";

        // shared independent variable
        fout << "independent_variables:\n";
        fout << "- header: {name: " << YAMLQuote(xName);
        if (!xUnits.empty()) fout << ", units: " << YAMLQuote(xUnits);
        fout << "}\n";
        fout << "  values:\n";

        const auto &refRows = groupedFormatted[g][0];
        for (const auto &row : refRows) {
            fout << "  - {high: " << row.xhigh << ", low: " << row.xlow << "}\n";
        }

        fout << "name: " << YAMLQuote(tableNamePrefix + " " + to_string(g + 1)) << "\n";
    }

    return true;
}

// ============================================================
// Main user function
// ============================================================

void FormatterTH1VtoHEPDataOverlay(
    const char *inFile,
    const char *outFile,
    const char *xName = "$p_{\\rm{T}}$",
    const char *xUnits = "GeV/c",
    const char *yUnits = "",
    const char *errorLabel = "total",
    const char *qualifierName = "",
    const char *qualifierValue = "",
    const char *descriptionPrefix = "Histograms for",
    const char *tableNamePrefix = "Table"
) {
    string file_name = inFile;
    if (!FileExists(file_name)) {
        cout << "Input file not found: " << file_name << "\n";
        return;
    }

    vector<TH1VBlock> blocks = ReadTH1VBlocks(file_name);
    if (blocks.empty()) {
        cout << "No valid TH1V blocks found.\n";
        return;
    }

    map<string, vector<size_t>> groups;
    for (size_t i = 0; i < blocks.size(); i++) {
        groups[MakeBinningKey(blocks[i])].push_back(i);
    }

    vector<vector<TH1VBlock>> groupedBlocks;
    vector<vector<vector<FormattedRow>>> groupedFormatted;

    for (const auto &kv : groups) {
        vector<TH1VBlock> gb;
        vector<vector<FormattedRow>> gf;

        for (size_t idx : kv.second) {
            gb.push_back(blocks[idx]);
            gf.push_back(FormatTH1VBlock(blocks[idx]));
        }

        groupedBlocks.push_back(gb);
        groupedFormatted.push_back(gf);
    }

    bool ok = WriteCombinedHEPDataDocuments(
        outFile,
        groupedBlocks,
        groupedFormatted,
        xName, xUnits, yUnits, errorLabel,
        qualifierName, qualifierValue,
        descriptionPrefix, tableNamePrefix
    );

    if (!ok) {
        cout << "Could not write output file: " << outFile << "\n";
        return;
    }

    cout << "Wrote combined HEPData YAML to: " << outFile << "\n";
    cout << "Number of overlaid tables: " << groupedBlocks.size() << "\n";
    for (size_t g = 0; g < groupedBlocks.size(); g++) {
        cout << "  Group " << g + 1 << ": ";
        for (size_t i = 0; i < groupedBlocks[g].size(); i++) {
            cout << groupedBlocks[g][i].hist_name;
            if (i + 1 < groupedBlocks[g].size()) cout << ", ";
        }
        cout << "\n";
    }
}

void FormatterTH1VtoHEPDataOverlayHelp() {
    cout << R"(
Usage:
  root -l
  .L formatter_TH1V_toHEPData_overlay.C+
  FormatterTH1VtoHEPDataOverlay("input.txt", "output.yaml");

Example:
  FormatterTH1VtoHEPDataOverlay(
    "input.txt",
    "output.yaml",
    "$p_{\rm{T}}$",
    "GeV/c",
    "",
    "total",
    "SQRT(S)/NUCLEON",
    "200 GeV",
    "Ratios for",
    "Table"
  );
)";
}
