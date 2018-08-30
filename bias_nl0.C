{

	const int n = 162;
	double pe[n],epe[n],tmp,e[n],ratio[n],sigma[n],esigma[n];
	double res[n],eres[n];
	double res2[n],bias[n],ebias[n];
	TF1 *fun0 = new TF1("fun0","(x*[0]+[1]*x*x+[2])/(1+[3]*exp(-[4]*x))",0,100);
	ifstream fin0("result");
//	ifstream fin0("/home/zhangfy/resolution/center_v3/result");

	double p0,p1,p2,p3,p4;
	fin0>>p0>>p1>>p2>>p3>>p4;
	fun0->SetParameters(p0,p1,p2,p3,p4);

	//TF1 *fun0 = new TF1("fun0","sqrt(2.993*2.993/x+0.7276*0.7276)",0.01,8);
	//TF1 *fun0 = new TF1("fun0","sqrt([0]*[0]/x+[1]*[1])",0.01,8);
	ifstream fin("positron_data");
	for(int i=0;i<n;i++){
		fin>>pe[i]>>epe[i]>>tmp>>sigma[i]>>esigma[i];
		e[i] = (i)/10.0/2;
		res[i] = 100*sigma[i]/pe[i];
		eres[i] = 100*esigma[i]/pe[i];
		double erec = pe[i]/1.3056438e3;
		double etrue = fun0->GetX(erec - 1259.53/1.3056438e3) + 1.022; 
		e[i] += 1.022;
		epe[i]/=1.3056438e3;
		bias[i] = 100*(etrue - e[i])/e[i];
		ebias[i] = 100*epe[i]/e[i];
		cout << etrue <<"\t"<< e[i] <<"\n";
	}
	TGraphErrors *t = new TGraphErrors(n,e,bias,0,ebias);
	t->GetXaxis()->SetRangeUser(1,8);
	t->Draw();
	
//	TGraphErrors *T = new TGraphErrors(n,pe,bias,0,ebias);
//	TGraphErrors *T = new TGraphErrors(n,pe,res,0,eres);
//		
//    	T->SetMarkerStyle(20);
//    	T->SetMarkerSize(1);
//    	T->Draw();
//    T->GetXaxis()->SetTitle("E_{rec} [MeV]");
//    T->GetYaxis()->SetTitle("resolution [%]");
//
//	fun0->Draw("same");
//
//    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
//    l->AddEntry(T,"data","p");
//    l->AddEntry(fun0,"prediction","l");
//    l->Draw();

}
