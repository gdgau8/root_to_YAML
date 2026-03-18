//
//  rootToText.c
//  
//
//  Created by Gabriel Dale-Gau on 05.01.2026.
//

#include <stdio.h>
#include <TFile.h>
#include <TKey.h>
#include <TObject.h>
#include <TClass.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include <TArrayD.h>

#include <fstream>
#include <iostream>
#include <string>

void rootToText(const char* infile  = "input.root",
                       const char* outfile = "output_varplus.txt")
{
    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening ROOT file: " << infile << "\n";
        return;
    }

    std::ofstream out(outfile);
    if (!out.is_open()) {
        std::cerr << "Error opening output text file: " << outfile << "\n";
        f->Close();
        return;
    }

    TIter next(f->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;

        // ---------------- TProfile (variable bins, ROOT 6.37-compatible) ----------------
        if (obj->InheritsFrom(TProfile::Class())) {
            TProfile* p = (TProfile*)obj;

            const TArrayD* bw2 = p->GetBinSumw2(); // sum of weights^2 per bin (if allocated)
            int hasW2 = (bw2 && bw2->GetSize() > 0) ? 1 : 0;

            const char* errOpt = p->GetErrorOption();
            if (!errOpt) errOpt = "";

            out << "# PROFILEV " << p->GetName() << " " << p->GetNbinsX()
                << " " << hasW2 << " " << errOpt << "\n";

            // Format per bin:
            // xLow xHigh sumY sumY2 entries sumW2
            for (int i = 1; i <= p->GetNbinsX(); ++i) {
                double lo = p->GetXaxis()->GetBinLowEdge(i);
                double hi = p->GetXaxis()->GetBinUpEdge(i);

                // Internals:
                // (*p)[i] is sum of Y in bin i
                // (*p->GetSumw2())[i] is sum of Y^2 in bin i (if present)
                double sumY  = (*p)[i];
                double sumY2 = 0.0;
                if (p->GetSumw2() && p->GetSumw2()->GetSize() > i) {
                    sumY2 = (*p->GetSumw2())[i];
                }

                // Entries is sum of weights / effective entries
                double entries = p->GetBinEntries(i);

                // Sum of weights^2, only if allocated
                double sumW2 = 0.0;
                if (hasW2 && bw2->GetSize() > i) {
                    sumW2 = (*bw2)[i];
                }

                out << lo << " " << hi << " "
                    << sumY << " " << sumY2 << " "
                    << entries << " " << sumW2 << "\n";
            }
            out << "\n";
            continue;
        }

        // ---------------- TH2 (variable bins) ----------------
        if (obj->InheritsFrom(TH2::Class())) {
            TH2* h = (TH2*)obj;

            out << "# TH2V " << h->GetName() << " "
                << h->GetNbinsX() << " " << h->GetNbinsY() << "\n";

            // Format per cell:
            // xLow xHigh yLow yHigh content error
            for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
                double xlo = h->GetXaxis()->GetBinLowEdge(ix);
                double xhi = h->GetXaxis()->GetBinUpEdge(ix);

                for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                    double ylo = h->GetYaxis()->GetBinLowEdge(iy);
                    double yhi = h->GetYaxis()->GetBinUpEdge(iy);

                    out << xlo << " " << xhi << " "
                        << ylo << " " << yhi << " "
                        << h->GetBinContent(ix, iy) << " "
                        << h->GetBinError(ix, iy) << "\n";
                }
            }
            out << "\n";
            continue;
        }

        // ---------------- TH1 (variable bins; exclude TH2 and TProfile) ----------------
        if (obj->InheritsFrom(TH1::Class()) &&
            !obj->InheritsFrom(TH2::Class()) &&
            !obj->InheritsFrom(TProfile::Class())) {

            TH1* h = (TH1*)obj;

            out << "# TH1V " << h->GetName() << " " << h->GetNbinsX() << "\n";
            // Format per bin:
            // xLow xHigh content error
            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                double lo = h->GetXaxis()->GetBinLowEdge(i);
                double hi = h->GetXaxis()->GetBinUpEdge(i);

                out << lo << " " << hi << " "
                    << h->GetBinContent(i) << " "
                    << h->GetBinError(i) << "\n";
            }
            out << "\n";
            continue;
        }

        // ---------------- Graphs ----------------
        if (obj->InheritsFrom(TGraphAsymmErrors::Class())) {
            TGraphAsymmErrors* g = (TGraphAsymmErrors*)obj;
            int n = g->GetN();
            out << "# GRAPHASYMM " << g->GetName() << " " << n << "\n";
            for (int i = 0; i < n; ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                out << x << " " << y << " "
                    << g->GetErrorXlow(i)  << " " << g->GetErrorXhigh(i) << " "
                    << g->GetErrorYlow(i)  << " " << g->GetErrorYhigh(i) << "\n";
            }
            out << "\n";
            continue;
        }

        if (obj->InheritsFrom(TGraphErrors::Class())) {
            TGraphErrors* g = (TGraphErrors*)obj;
            int n = g->GetN();
            out << "# GRAPHERR " << g->GetName() << " " << n << "\n";
            for (int i = 0; i < n; ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                out << x << " " << y << " "
                    << g->GetErrorX(i) << " " << g->GetErrorY(i) << "\n";
            }
            out << "\n";
            continue;
        }

        if (obj->InheritsFrom(TGraph::Class())) {
            TGraph* g = (TGraph*)obj;
            int n = g->GetN();
            out << "# GRAPH " << g->GetName() << " " << n << "\n";
            for (int i = 0; i < n; ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                out << x << " " << y << "\n";
            }
            out << "\n";
            continue;
        }

        // else ignore other object types
    }

    out.close();
    f->Close();
    std::cout << "Wrote objects to " << outfile << "\n";
}
