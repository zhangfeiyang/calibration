{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Aug 28 16:07:30 2018) by ROOT version5.34/11
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   c1->SetHighLightColor(2);
   c1->Range(-0.8364876,0.9023977,7.528388,1.101891);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(8);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(0.8);
   gre->SetPoint(0,0.5109989,0.9390517);
   gre->SetPointError(0,0,0.003405171);
   gre->SetPoint(1,0.6617,0.9545073);
   gre->SetPointError(1,0,0.003109317);
   gre->SetPoint(2,0.83485,0.9805087);
   gre->SetPointError(2,0,0.003290862);
   gre->SetPoint(3,1.252868,1.000109);
   gre->SetPointError(3,0,0.003247504);
   gre->SetPoint(4,1.4608,1.005788);
   gre->SetPointError(4,0,0.00311754);
   gre->SetPoint(5,2.223,1.031773);
   gre->SetPointError(5,0,0.003099322);
   gre->SetPoint(6,4.95,1.048101);
   gre->SetPointError(6,0,0.003161095);
   gre->SetPoint(7,6.13,1.063898);
   gre->SetPointError(7,0,0.004743921);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,6.6919);
   Graph_Graph1->SetMinimum(0.922347);
   Graph_Graph1->SetMaximum(1.081941);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   
   TF1 *f1 = new TF1("*f1",0,8,5);
    //The original function :  had originally been created by:
    //TF1 *f1 = new TF1("f1",,0,8,5);
   f1->SetRange(0,8);
   f1->SetName("f1");
   f1->SetTitle("");
   f1->SetSavedPoint(0,0.9379412);
   f1->SetSavedPoint(1,0.9379412);
   f1->SetSavedPoint(2,0.9379412);
   f1->SetSavedPoint(3,0.9379412);
   f1->SetSavedPoint(4,0.9379412);
   f1->SetSavedPoint(5,0.9379412);
   f1->SetSavedPoint(6,0.9379412);
   f1->SetSavedPoint(7,0.9379412);
   f1->SetSavedPoint(8,0.9585321);
   f1->SetSavedPoint(9,0.9585321);
   f1->SetSavedPoint(10,0.9750614);
   f1->SetSavedPoint(11,0.9750614);
   f1->SetSavedPoint(12,0.9750614);
   f1->SetSavedPoint(13,0.9750614);
   f1->SetSavedPoint(14,1.00047);
   f1->SetSavedPoint(15,1.00047);
   f1->SetSavedPoint(16,1.00047);
   f1->SetSavedPoint(17,1.008939);
   f1->SetSavedPoint(18,1.008939);
   f1->SetSavedPoint(19,1.008939);
   f1->SetSavedPoint(20,1.008939);
   f1->SetSavedPoint(21,1.008939);
   f1->SetSavedPoint(22,1.008939);
   f1->SetSavedPoint(23,1.008939);
   f1->SetSavedPoint(24,1.029202);
   f1->SetSavedPoint(25,1.029202);
   f1->SetSavedPoint(26,1.029202);
   f1->SetSavedPoint(27,1.029202);
   f1->SetSavedPoint(28,1.029202);
   f1->SetSavedPoint(29,1.029202);
   f1->SetSavedPoint(30,1.029202);
   f1->SetSavedPoint(31,1.029202);
   f1->SetSavedPoint(32,1.029202);
   f1->SetSavedPoint(33,1.029202);
   f1->SetSavedPoint(34,1.029202);
   f1->SetSavedPoint(35,1.029202);
   f1->SetSavedPoint(36,1.029202);
   f1->SetSavedPoint(37,1.029202);
   f1->SetSavedPoint(38,1.029202);
   f1->SetSavedPoint(39,1.029202);
   f1->SetSavedPoint(40,1.029202);
   f1->SetSavedPoint(41,1.029202);
   f1->SetSavedPoint(42,1.029202);
   f1->SetSavedPoint(43,1.029202);
   f1->SetSavedPoint(44,1.029202);
   f1->SetSavedPoint(45,1.049708);
   f1->SetSavedPoint(46,1.049708);
   f1->SetSavedPoint(47,1.049708);
   f1->SetSavedPoint(48,1.049708);
   f1->SetSavedPoint(49,1.049708);
   f1->SetSavedPoint(50,1.049708);
   f1->SetSavedPoint(51,1.049708);
   f1->SetSavedPoint(52,1.049708);
   f1->SetSavedPoint(53,1.049708);
   f1->SetSavedPoint(54,1.049708);
   f1->SetSavedPoint(55,1.049708);
   f1->SetSavedPoint(56,1.049708);
   f1->SetSavedPoint(57,1.049708);
   f1->SetSavedPoint(58,1.049708);
   f1->SetSavedPoint(59,1.049708);
   f1->SetSavedPoint(60,1.049708);
   f1->SetSavedPoint(61,1.049708);
   f1->SetSavedPoint(62,1.049708);
   f1->SetSavedPoint(63,1.049708);
   f1->SetSavedPoint(64,1.049708);
   f1->SetSavedPoint(65,1.049708);
   f1->SetSavedPoint(66,1.049708);
   f1->SetSavedPoint(67,1.049708);
   f1->SetSavedPoint(68,1.049708);
   f1->SetSavedPoint(69,1.049708);
   f1->SetSavedPoint(70,1.06234);
   f1->SetSavedPoint(71,1.06234);
   f1->SetSavedPoint(72,1.06234);
   f1->SetSavedPoint(73,1.06234);
   f1->SetSavedPoint(74,1.06234);
   f1->SetSavedPoint(75,1.06234);
   f1->SetSavedPoint(76,1.06234);
   f1->SetSavedPoint(77,1.06234);
   f1->SetSavedPoint(78,1.06234);
   f1->SetSavedPoint(79,1.06234);
   f1->SetSavedPoint(80,1.06234);
   f1->SetSavedPoint(81,1.06234);
   f1->SetSavedPoint(82,1.06234);
   f1->SetSavedPoint(83,1.06234);
   f1->SetSavedPoint(84,1.06234);
   f1->SetSavedPoint(85,1.06234);
   f1->SetSavedPoint(86,1.067876);
   f1->SetSavedPoint(87,1.067876);
   f1->SetSavedPoint(88,1.067876);
   f1->SetSavedPoint(89,1.067876);
   f1->SetSavedPoint(90,1.067876);
   f1->SetSavedPoint(91,1.067876);
   f1->SetSavedPoint(92,1.067876);
   f1->SetSavedPoint(93,1.067876);
   f1->SetSavedPoint(94,1.067876);
   f1->SetSavedPoint(95,1.067876);
   f1->SetSavedPoint(96,1.067876);
   f1->SetSavedPoint(97,1.067876);
   f1->SetSavedPoint(98,1.067876);
   f1->SetSavedPoint(99,1.067876);
   f1->SetSavedPoint(100,1.067876);
   f1->SetSavedPoint(101,0);
   f1->SetSavedPoint(102,8);
   f1->SetFillColor(19);
   f1->SetFillStyle(0);
   f1->SetLineColor(2);
   f1->SetLineWidth(2);
   f1->SetChisquare(6.610111);
   f1->SetNDF(3);
   f1->GetXaxis()->SetLabelFont(42);
   f1->GetXaxis()->SetLabelSize(0.035);
   f1->GetXaxis()->SetTitleSize(0.035);
   f1->GetXaxis()->SetTitleFont(42);
   f1->GetYaxis()->SetLabelFont(42);
   f1->GetYaxis()->SetLabelSize(0.035);
   f1->GetYaxis()->SetTitleSize(0.035);
   f1->GetYaxis()->SetTitleFont(42);
   f1->SetParameter(0,1.064513);
   f1->SetParError(0,0.02918158);
   f1->SetParLimits(0,0,0);
   f1->SetParameter(1,0.00336619);
   f1->SetParError(1,0.006798957);
   f1->SetParLimits(1,0,0);
   f1->SetParameter(2,-0.002462761);
   f1->SetParError(2,0.001030202);
   f1->SetParLimits(2,0,0);
   f1->SetParameter(3,0.0578301);
   f1->SetParError(3,0.04229992);
   f1->SetParLimits(3,0,0);
   f1->SetParameter(4,1.440251);
   f1->SetParError(4,1.806371);
   f1->SetParLimits(4,0,0);
   gre->GetListOfFunctions()->Add(f1);
   gre->Draw("ap");
   
   TGraph *graph = new TGraph(8);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(0.8);
   graph->SetPoint(0,0.5109989,0.9379412065);
   graph->SetPoint(1,0.6617,0.9585320835);
   graph->SetPoint(2,0.83485,0.9750613515);
   graph->SetPoint(3,1.2528685,1.000470472);
   graph->SetPoint(4,1.4608,1.008939179);
   graph->SetPoint(5,2.223,1.029202377);
   graph->SetPoint(6,4.95,1.049707591);
   graph->SetPoint(7,6.13,1.062340339);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,6.6919);
   Graph_Graph1->SetMinimum(0.9255013);
   Graph_Graph1->SetMaximum(1.07478);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("pc");
   
   TPaveText *pt = new TPaveText(0.4397126,0.9339831,0.5602874,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("Graph");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
