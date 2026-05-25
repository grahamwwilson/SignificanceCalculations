{

// For the range 0.05 to 5.2

   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",1200,1000);

   c1->SetFillColor(0);
   c1->SetGrid();
   c1->SetTicks(1,1);
   c1->GetFrame()->SetFillColor(0);
   c1->GetFrame()->SetBorderSize(12);
   c1->SetLeftMargin(0.14);
   c1->SetRightMargin(0.08);

   TLatex *title = new TLatex(500.0,1.8e4,"Bhabha and Dimuon Cross-sections");
   title->SetTextSize(0.05);
   title->SetTextFont(62);
   title->SetTextAlign(22);   // centered horizontally and vertically

   Size_t msize=4.0;
   int width=5;
   int dw=2;

   // By default expects x,y,dx,dy
   TGraph *gr = new TGraph("ZScore-Calculator-SD-PAll.tpl","%lg %lg");
   gr->SetTitle("Background dependence of std. dev. (numerical estimate)");
//   gr->SetName("0.125 mrad");
   gr->SetMarkerSize(0.20*msize);
   gr->SetMarkerColor(4);
   gr->SetLineColor(7);
   gr->SetLineWidth(width-dw);
   gr->SetMarkerStyle(20);
   gr->GetXaxis()->SetTitle("#hat{#mu}_{b}");
   gr->GetYaxis()->SetTitle("Std. dev. of #it{Z}^{S}_{N}");
//   TAxis xaxis=gr->GetXaxis();
   gr->GetXaxis()->SetLimits(0.05, 5.2);
   gr->SetMaximum(1.1);
   gr->SetMinimum(0.40);
   
   c1->SetLogx();

/*
   gr3 = new TGraphErrors("sense_all_PLR.dat"); 
   gr3->SetMarkerSize(msize);
   gr3->SetMarkerColor(2);
   gr3->SetLineColor(2);
   gr3->SetLineWidth(width-dw);   
   gr3->SetMarkerStyle(22);

   gr4 = new TGraphErrors("sense_all_LR.dat"); 
   gr4->SetMarkerSize(msize);
   gr4->SetMarkerColor(4);
   gr4->SetLineColor(4);
   gr4->SetLineWidth(width-dw);   
   gr4->SetMarkerStyle(20);
*/   

   // make a TMultiGraph (needed for legend management) and draw it
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gr,"apc");   
//   mg->Add(gr4,"cp");   
//   mg->Add(gr3,"cp");  

   mg->Draw("");

   title->Draw();

//   TLegend* leg = new TLegend(0.65, 0.55, 0.9, 0.95,"Angular Resolution","brNDC");
//   TLegend* leg = new TLegend(0.30, 0.70, 0.9, 0.88);

//   leg->AddEntry(gr, "480 1/10-X_{0} 750 um Si, 500 um G10, W ", "P");
//   leg->AddEntry(gr3, "P(e-, e+) = (-0.8, +0.3)", "P");
//   leg->AddEntry(gr4, "P(e-, e+) = (-1, +1)", "P");         
   

//   leg->SetBorderSize(0);
//   leg->Draw();

   TLatex *t = new TLatex();
   t->SetTextFont(42);
   t->SetTextColor(kBlack);
   t->SetTextSize(0.04);
   t->SetTextAlign(12);
//   t->DrawLatex(0.5,3.5,"#pm 1 #sigma band (#chi^{2}/#nu = 2.3/4)");
   
   TLatex *t2 = new TLatex();
   t2->SetTextFont(42);
   t2->SetTextColor(kMagenta);
   t2->SetTextSize(0.04);
   t2->SetTextAlign(12);
   t2->DrawLatex(0.13,1.05,"Evaluated in the Poisson limit (#hat{#sigma}_{b} = 0)");   

   c1->Update();
   
   double y0 = 3.66364;
   double dy0 = 0.007673;
   double yminus = y0 - dy0;
   double yplus  = y0 + dy0;
   TLine* line1 = new TLine(-0.2, yminus, 2.6, yminus);
   TLine* line2 = new TLine(-0.2, yplus, 2.6, yplus); 
   TLine* line0 = new TLine(-0.2, 0.5*(yplus+yminus), 2.6, 0.5*(yplus+yminus));
   line1->SetLineColor(8);
   line1->SetLineWidth(3);
   line2->SetLineColor(8);
   line2->SetLineWidth(3);        
   line0->SetLineWidth(3);
//   line1->Draw();
//   line2->Draw();  
//   line0->Draw();
}
