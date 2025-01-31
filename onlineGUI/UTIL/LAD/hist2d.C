void hist2d(TString hist1name, TString hist2name, TString label1 = "", TString label2 = "", TString htitle = "", Bool_t golden = false) {
    // Retrieve histograms from the current directory
    TH1F* H1 = (TH1F*)gDirectory->Get(hist1name);
    TH1F* H2 = (TH1F*)gDirectory->Get(hist2name);

    if (H1 && H2) {
        // Get number of bins from histograms (assuming they have the same binning)
        Int_t nbinsX = H1->GetNbinsX();

        // Create a TH2F to plot H1 vs H2
        TH2F* H2D = new TH2F("H1_vs_H2", htitle, nbinsX, H1->GetXaxis()->GetXmin(), H1->GetXaxis()->GetXmax(),
                                            nbinsX, H2->GetXaxis()->GetXmin(), H2->GetXaxis()->GetXmax());

        // Set axis titles
        H2D->GetXaxis()->SetTitle(label1);
        H2D->GetYaxis()->SetTitle(label2);

        // Fill TH2F with the (x, y) pairs from H1 and H2 bin contents
        for (int i = 1; i <= nbinsX; i++) {
            double x_val = H1->GetBinContent(i);
            double y_val = H2->GetBinContent(i);
            H2D->Fill(H1->GetXaxis()->GetBinCenter(i), y_val);
        }

        // Apply style options based on "golden" flag
        if (golden) {
            H2D->SetMarkerColor(30);
            H2D->SetMarkerStyle(20);
            H2D->SetMarkerSize(0.8);
        } else {
            H2D->SetMarkerColor(kBlue);
            H2D->SetMarkerStyle(21);
            H2D->SetMarkerSize(0.8);
        }

        // Draw the 2D histogram
        H2D->Draw("COLZ");

    } else {
        if (!H1) cout << "Histogram " << hist1name << " does not exist" << endl;
        if (!H2) cout << "Histogram " << hist2name << " does not exist" << endl;
    }
}
