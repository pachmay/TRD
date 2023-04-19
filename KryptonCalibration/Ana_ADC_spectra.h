
using namespace std;

#include "TString.h"
#include "TH2F.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveText.h"
//#include <iostream.h>
#include <fstream>
#include <sstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"
#include "TList.h"
#include "TChain.h"

#include "TMinuit.h"
#include "TFitter.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"

#include "Minuit2/Minuit2Minimizer.h"

//#include "TPython.h"
#include "TArrow.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TNtuple.h"
#include "TDatime.h"
#include "TClonesArray.h"

//#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <string>     // std::string, std::stoi

#include "bad_chambers.h"


//------------------------------------------------------------------------------------------------------------
static const Float_t Pi = TMath::Pi();
static TRandom ran;
static TString HistName, HistNameB;
static char NoP[50];
static TF1* func_Gauss_fit;
static TF1* func_Poly_fit;
static TF1* func_Krypton_fit;
static TF1* func_Krypton_Student_t_fit;
static TF1* func_Landau_fit;
static TF1* func_Expo_fit;
static TF1* func_Student_t_fit;
static const Int_t N_TRD = 540;
static const Int_t N_TRD_sectors = 18;
static const Int_t N_TRD_layers  = 6;
static const Int_t N_TRD_stacks  = 5;
static const Int_t N_TRD_rows    = 16;
static const Int_t N_TRD_columns = 144;
static Int_t merge_N_rows    = 2;
static Int_t merge_N_columns = 16;
static TGraph* tg_pressure_vs_time;
static TGraph* tg_pressure_vs_run_id;

static Int_t N_columns_merge;
static Int_t N_rows_merge;
static Int_t N_run_ids;
static Int_t det_full;
static Int_t sector_plot;
static vector< vector< vector<TH1F*> > > vec_ADC_pads;
static vector< vector< vector< vector<TH1F*> > > > vec_ADC_pads_corr;
static vector< vector< vector<TH1F*> > > vec_ADC_pads_merge;
static vector< vector< vector< vector<TH1F*> > > > vec_ADC_pads_merge_corr;
static vector< vector<TH1F*> > vec_ADC_pads_merge_sum_run_ids;
static vector< vector<TH1F*> > vec_ADC_pads_merge_sum_run_ids_corr;
static vector< vector< vector< vector<Double_t> > > > vec_merge_fit_par;
static vector< vector< vector< vector<Double_t> > > > vec_merge_mean_pressure_fit_par;
static vector< vector< vector<TGraphErrors*> > > tg_mean_vs_pressure;
static vector<TH1F*> h_slopes;
static vector<TH1F*> h_slopes_over_b;
static vector<Double_t> vec_pressure;
static vector<TH1F*> vec_h_ADC_TRD_chambers;
static Int_t arr_color_row_merge[8] = {kBlack,kRed,kBlue,kGreen+1,kMagenta,kCyan,kOrange+1,kTeal+8};
static Double_t pressure_ref = 965.0; // reference pressure, used for correcting the ADC pressure dependence
static TH1F* h_corr_ADC_full;
static vector<TCanvas*> vec_can_TRD_det_ADC_sectors;
static TFile* inputfile;
static vector<UInt_t> vec_run_ids;
static TH1F* h_ADC_vs_pressure_fit_par0;
static TH1F* h_ADC_vs_pressure_fit_par1;
static TH2F* h2D_ADC_vs_merged_pads;
static TH2F* h2D_ADC_vs_merged_pads_all;
static TCanvas* can_ADC_fit;
static TString main_data_dir;
static vector<TCanvas*> vec_c_slopes;
static vector<TCanvas*> vec_c_mean_vs_pressure;
static vector<TH1F*> vec_h_peak_params;
static TH2F* h2D_ADC_mean_main_Krypton_xy_TRD;
static TH2F* h2D_ADC_sigma_main_Krypton_xy_TRD;
static TH2F* h2D_ADC_amplitude_main_Krypton_xy_TRD;
static TH2F* h2D_ADC_chi2_main_Krypton_xy_TRD;
static TH2F* h2D_ADC_mean_second_Krypton_xy_TRD;
static TH2F* h2D_ADC_sigma_second_Krypton_xy_TRD;
static TH2F* h2D_ADC_amplitude_second_Krypton_xy_TRD;
static TH2F* h2D_ADC_chi2_second_Krypton_xy_TRD;
static TCanvas* can_ADC;
static TFile* file_Krypton;
static TNtuple* nt_Krypton;
static TH1F* h_ADC_full_merge;
static vector<TH1F*> h_ADC_full_merge_det;
static vector<Int_t> vec_bad_chambers     = vec_bad_chambers_init;
static vector<Int_t> vec_add_bad_chambers = vec_add_bad_chambers_init;

vector< TProfile* > vec_h_fit_params;
vector< TProfile* > vec_h_fit_params_merged;
TCanvas* can_ADC_fit_pretty;
//------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,
                      Float_t size=0.06,Int_t color=1,Float_t angle=0.0,
                      Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color
    // align: 1 left aligned, 32, right aligned

    // align = 10*HorizontalAlign + VerticalAlign
    // For horizontal alignment the following convention applies:
    // 1=left adjusted, 2=centered, 3=right adjusted
    // For vertical alignment the following convention applies:
    // 1=bottom adjusted, 2=centered, 3=top adjusted

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Krypton_fit(Double_t* x_val, Double_t* par)
{
    // With pedestal
    Double_t x, y, par0, par1, par2, par3, par4;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    par3  = par[3];
    par4  = par[4];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0)
        + par0*(0.155/0.58)*TMath::Student((x-par1*(29.0/41.6))/(par2*(42.0/41.6)),par3)
        + par0*(0.04/0.58)*TMath::Student((x-par1*(19.6/41.6))/(par2*(19.6/41.6)),par3)
        + par0*(0.25/0.58)*TMath::Student((x-par1*(12.6/41.6))/(par2*(35.0/41.6)),par3)
        + par0*(0.08/0.58)*TMath::Student((x-par1*(9.6/41.6))/(par2*(9.6/41.6)),par3)
        + par3 + par4*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Krypton_Student_t_fit(Double_t* x_val, Double_t* par)
{
    // With pedestal
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0*TMath::Student((x-par1)/par2,par3)
        //+ par0*(0.14/0.58)*TMath::Student((x-par1*(29.0/41.6))/(par2*(29.0/41.6)),par3)
        //+ par0*(0.04/0.58)*TMath::Student((x-par1*(19.6/41.6))/(par2*(19.6/41.6)),par3)
        //+ par0*(0.15/0.58)*TMath::Student((x-par1*(12.6/41.6))/(par2*(12.6/41.6)),par3)
        //+ par0*(0.12/0.58)*TMath::Student((x-par1*(9.6/41.6))/(par2*(9.6/41.6)),par3)

        + par0*(0.155/0.58)*TMath::Student((x-par1*(29.0/41.6))/(par2*(42.0/41.6)),par3)
        + par0*(0.04/0.58)*TMath::Student((x-par1*(19.6/41.6))/(par2*(19.6/41.6)),par3)
        + par0*(0.25/0.58)*TMath::Student((x-par1*(12.6/41.6))/(par2*(35.0/41.6)),par3)
        + par0*(0.08/0.58)*TMath::Student((x-par1*(9.6/41.6))/(par2*(9.6/41.6)),par3)

        + par4 + par5*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LandauFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0]; // amplitude
    par1  = par[1]; // mean
    par2  = par[2]; // sigma
    x = x_val[0];
    y = par0*TMath::Landau(x,par1,par2);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t ExpoFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, Amplitude, T, mu;
    Amplitude  = par[0];
    T          = par[1];
    mu         = par[2];
    x          = x_val[0];
    y          = Amplitude*TMath::Exp(-(x-mu)/T);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Student_t_FitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3;
    par0  = par[0]; // amplitude
    par1  = par[1]; // mean
    par2  = par[2]; // width
    par3  = par[3]; // NDF
    x = x_val[0];
    y = par0*TMath::Student((x-par1)/par2,par3);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    par3  = par[3];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + par3;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Init_functions()
{
    cout << "Fit functions initialized" << endl;
    func_Gauss_fit             = new TF1("func_Gauss_fit",GaussFitFunc,0,15000,4);
    func_Gauss_fit->SetNpx(400);
    func_Poly_fit              = new TF1("func_Poly_fit",PolyFitFunc,0,15000,6);
    func_Krypton_fit           = new TF1("func_Krypton_fit",Krypton_fit,0,15000,5);
    func_Landau_fit            = new TF1("func_Landau_fit",LandauFitFunc,0,15000,3);
    func_Expo_fit              = new TF1("func_Expo_fit",ExpoFitFunc,0,15000,2);
    func_Student_t_fit         = new TF1("func_Student_t_fit",Student_t_FitFunc,0,15000,4);
    func_Krypton_Student_t_fit = new TF1("func_Krypton_Student_t_fit",Krypton_Student_t_fit,0,15000,6);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
class Class_TRD_Krypton_ADC : public TObject
{
private:
    //
    UShort_t   ADC;
    UInt_t     channel;
    UShort_t   rmsTime;

public:
    Class_TRD_Krypton_ADC() :
        ADC(0), channel(0), rmsTime(0)
    {
    }
    ~Class_TRD_Krypton_ADC()
    {
    }

    // setters
    void set_ADC(UShort_t us)                             {ADC               = us;}
    void set_channel(UInt_t i)                            {channel           = i;}
    void set_rmsTime(UShort_t us)                         {rmsTime           = us;}


    // getters
    UShort_t  get_ADC() const                             {return ADC;}
    UInt_t    get_channel() const                         {return channel;}
    Double_t  get_rmsTime() const                         {return ((Double_t)rmsTime)/10000.0;}


    ClassDef(Class_TRD_Krypton_ADC,1);
};




class Class_TRD_Krypton : public TObject
{
private:

    UInt_t        run_id;
    UShort_t      fNumADCs;

    // The-> behind TClonesArray is important!
    TClonesArray* fTRDADC;      //->

public:
    Class_TRD_Krypton() :
        run_id(0),fNumADCs(0),fTRDADC(0)
    {
        fTRDADC = new TClonesArray( "Class_TRD_Krypton_ADC", 10 );
    }
        ~Class_TRD_Krypton()
        {
            delete fTRDADC;
            fTRDADC = NULL;
        }

        // setters
        void      set_run_id(UInt_t i)         {run_id = i; }

        // getters
        UInt_t    get_run_id() const           {return run_id;}

        Class_TRD_Krypton_ADC* createADC()
        {
            if (fNumADCs == fTRDADC->GetSize())
                fTRDADC->Expand( fNumADCs + 10 );
            if (fNumADCs >= 60000)
            {
                Fatal( "Class_TRD_Krypton::createADC()", "ERROR: Too many ADCs (>10000)!" );
                exit( 2 );
            }

            new((*fTRDADC)[fNumADCs++]) Class_TRD_Krypton_ADC;
            return (Class_TRD_Krypton_ADC*)((*fTRDADC)[fNumADCs - 1]);
        }
        void clearADCList()
        {
            fNumADCs   = 0;
            fTRDADC   ->Clear();
        }
        UShort_t getNumADCs() const
        {
            return fNumADCs;
        }
        Class_TRD_Krypton_ADC* getADC(UShort_t i) const
        {
            return i < fNumADCs ? (Class_TRD_Krypton_ADC*)((*fTRDADC)[i]) : NULL;
        }

        ClassDef(Class_TRD_Krypton,1);
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_2D_histo_and_canvas(TH2F* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.22);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    if(max_val > min_val)
    {
        hist->GetZaxis()->SetRangeUser(min_val,max_val);
    }
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Set_map_style_and_ranges(TH2F* h2D_map_in,TString x_title, TString y_title, TString z_title,
                              Double_t x_start_plot, Double_t x_stop_plot,
                              Double_t y_start_plot, Double_t y_stop_plot,
                              Double_t min_z_fraction_cut, Double_t max_z_fraction_cut
                             )
{
    h2D_map_in->SetStats(0);
    h2D_map_in->SetTitle("");
    h2D_map_in->GetXaxis()->SetTitleOffset(1.1);
    h2D_map_in->GetYaxis()->SetTitleOffset(1.4);
    h2D_map_in->GetZaxis()->SetTitleOffset(1.4);
    h2D_map_in->GetXaxis()->SetLabelSize(0.055);
    h2D_map_in->GetYaxis()->SetLabelSize(0.055);
    h2D_map_in->GetZaxis()->SetLabelSize(0.055);
    h2D_map_in->GetXaxis()->SetTitleSize(0.055);
    h2D_map_in->GetYaxis()->SetTitleSize(0.055);
    h2D_map_in->GetZaxis()->SetTitleSize(0.055);
    h2D_map_in->GetXaxis()->SetNdivisions(505,'N');
    h2D_map_in->GetYaxis()->SetNdivisions(505,'N');
    h2D_map_in->GetZaxis()->SetNdivisions(505,'N');
    h2D_map_in->GetXaxis()->CenterTitle();
    h2D_map_in->GetYaxis()->CenterTitle();
    h2D_map_in->GetZaxis()->CenterTitle();
    h2D_map_in->GetXaxis()->SetTitle(x_title.Data());
    h2D_map_in->GetYaxis()->SetTitle(y_title.Data());
    h2D_map_in->GetZaxis()->SetTitle(z_title.Data());

    Double_t xy_plot_range[2][2] =
    {
        {x_start_plot,y_start_plot}, // start, x,y
        {x_stop_plot,y_stop_plot} // stop, x, y
    };
    Int_t xy_bin_range[2][2];
    xy_bin_range[0][0] = h2D_map_in->GetXaxis()->FindBin(xy_plot_range[0][0]);
    xy_bin_range[1][0] = h2D_map_in->GetXaxis()->FindBin(xy_plot_range[1][0]);
    xy_bin_range[0][1] = h2D_map_in->GetYaxis()->FindBin(xy_plot_range[0][1]);
    xy_bin_range[1][1] = h2D_map_in->GetYaxis()->FindBin(xy_plot_range[1][1]);

    Double_t min_max_val[2] = {1000000.0,-1000000.0};
    for(Int_t i_bin_x = xy_bin_range[0][0]; i_bin_x <= xy_bin_range[1][0]; i_bin_x++)
    {
        for(Int_t i_bin_y = xy_bin_range[0][1]; i_bin_y <= xy_bin_range[1][1]; i_bin_y++)
        {
            Double_t bin_cont = h2D_map_in->GetBinContent(i_bin_x,i_bin_y);
            if(bin_cont < min_max_val[0]) min_max_val[0] = bin_cont;
            if(bin_cont > min_max_val[1]) min_max_val[1] = bin_cont;
        }
    }

    printf("min max vals: {%f,%f} \n",min_max_val[0],min_max_val[1]);

    TH1F* h_height = new TH1F("h_height","h_height",500,min_max_val[0],min_max_val[1]);
    for(Int_t i_bin_x = xy_bin_range[0][0]; i_bin_x <= xy_bin_range[1][0]; i_bin_x++)
    {
        for(Int_t i_bin_y = xy_bin_range[0][1]; i_bin_y <= xy_bin_range[1][1]; i_bin_y++)
        {
            Double_t bin_cont = h2D_map_in->GetBinContent(i_bin_x,i_bin_y);
            h_height->Fill(bin_cont);
        }
    }
    Double_t int_h_height = h_height->Integral(1,-1);
    if(int_h_height > 0.0)
    {
        Double_t bin_sum = 0.0;
        for(Int_t i_bin_x = 1; i_bin_x < h_height->GetNbinsX(); i_bin_x++)
        {
            Double_t bin_cont = h_height->GetBinContent(i_bin_x);
            Double_t bin_cent = h_height->GetBinCenter(i_bin_x);
            bin_sum += bin_cont;
            if((bin_sum/int_h_height) < min_z_fraction_cut)
            {
                min_max_val[0] = bin_cent;
            }
            if((bin_sum/int_h_height) > max_z_fraction_cut)
            {
                min_max_val[1] = bin_cent;
            }
        }
    }

    printf("min max vals cut: {%f,%f} \n",min_max_val[0],min_max_val[1]);

    h_height->Delete();

    h2D_map_in->GetXaxis()->SetRangeUser(xy_plot_range[0][0],xy_plot_range[1][0]);
    h2D_map_in->GetYaxis()->SetRangeUser(xy_plot_range[0][1],xy_plot_range[1][1]);
    h2D_map_in->GetZaxis()->SetRangeUser(min_max_val[0],min_max_val[1]);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_graph_and_canvas(TGraph* graph, TString x_title, TString y_title, Int_t color, Int_t style, Double_t size,
                               TString name, Int_t x_size, Int_t y_size,
                               Double_t x_min, Double_t x_max, Double_t y_min, Double_t y_max)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.1);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    TString hist_name = name;
    hist_name += "dummy";
    TH1F* hist = new TH1F(hist_name.Data(),hist_name.Data(),1000,x_min,x_max);
    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->SetTitle(x_title.Data());
    hist->GetYaxis()->SetTitle(y_title.Data());
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(y_min != y_max) hist->GetYaxis()->SetRangeUser(y_min,y_max);
    hist->SetLineColor(10);
    hist->DrawCopy("h");

    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(style);
    graph->SetMarkerSize(size);
    graph->Draw("p");

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Plot_pressure_vs_run_id()
{
    Double_t x_start = 0.0;
    Double_t x_stop  = 1.0;
    Double_t y_start = 0.0;
    Double_t y_stop  = 1.0;
    for(Int_t i_point = 0; i_point < (Int_t)tg_pressure_vs_run_id->GetN(); i_point++)
    {
        Double_t run_id, pressure;
        tg_pressure_vs_run_id->GetPoint(i_point,run_id,pressure);
        if(i_point == 0)
        {
            x_start = run_id;
            x_stop  = run_id;
            y_start = pressure;
            y_stop  = pressure;
        }
        else
        {
            if(run_id < x_start) x_start = run_id;
            if(run_id > x_stop)  x_stop  = run_id;
            if(pressure < y_start) y_start = pressure;
            if(pressure > y_stop)  y_stop  = pressure;
        }
        //printf("i_point: %d, run_id: %f, pressure: %f \n",i_point,run_id,pressure);
    }
    Draw_graph_and_canvas(tg_pressure_vs_run_id,"run id","pressure (hPa)",kBlack,20,0.6,"pressure_vs_run_id",900,500,x_start,x_stop,y_start,y_stop);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void init_ADC_pad_merge(Int_t sector_use)
{
    cout << "" << endl;
    printf("init_ADC_pad_merge for sector: %d \n",sector_use);
    N_run_ids = 1;
    vec_run_ids.push_back(0);

    nt_Krypton = new TNtuple("nt_Krypton","nt_Krypton","det:col:row:chi2:mean:sigma:amplitude");

    h_ADC_vs_pressure_fit_par0 = new TH1F("h_ADC_vs_pressure_fit_par0","h_ADC_vs_pressure_fit_par0",N_TRD*N_rows_merge*N_columns_merge,0,N_TRD*N_rows_merge*N_columns_merge);
    h_ADC_vs_pressure_fit_par1 = new TH1F("h_ADC_vs_pressure_fit_par1","h_ADC_vs_pressure_fit_par1",N_TRD*N_rows_merge*N_columns_merge,0,N_TRD*N_rows_merge*N_columns_merge);

    h2D_ADC_vs_merged_pads = new TH2F("h2D_ADC_vs_merged_pads","h2D_ADC_vs_merged_pads",N_columns_merge*(N_TRD_sectors-1)+N_columns_merge,0,N_columns_merge*(N_TRD_sectors-1)+N_columns_merge,(N_TRD_stacks-1)*N_TRD_layers*N_rows_merge+(N_TRD_layers-1)*N_rows_merge+N_rows_merge,0,(N_TRD_stacks-1)*N_TRD_layers*N_rows_merge+(N_TRD_layers-1)*N_rows_merge+N_rows_merge);

    h2D_ADC_vs_merged_pads_all= new TH2F("h2D_ADC_vs_merged_pads_all","h2D_ADC_vs_merged_pads_all",N_columns_merge*(N_TRD_sectors-1)+N_columns_merge,0,N_columns_merge*(N_TRD_sectors-1)+N_columns_merge,(N_TRD_stacks-1)*N_TRD_layers*N_rows_merge+(N_TRD_layers-1)*N_rows_merge+N_rows_merge,0,(N_TRD_stacks-1)*N_TRD_layers*N_rows_merge+(N_TRD_layers-1)*N_rows_merge+N_rows_merge);


    can_ADC_fit = new TCanvas("can_ADC_fit","can_ADC_fit",10,10,1400,1000);
    can_ADC_fit->SetFillColor(10);
    can_ADC_fit->SetTopMargin(0.05);
    can_ADC_fit->SetBottomMargin(0.2);
    can_ADC_fit->SetRightMargin(0.05);
    can_ADC_fit->SetLeftMargin(0.2);
    can_ADC_fit->SetTicks(1,1);
    can_ADC_fit->SetGrid(0,0);
    can_ADC_fit->Divide(N_columns_merge,N_rows_merge);

    can_ADC_fit_pretty = new TCanvas( "can_ADC_fit_pretty", "can_ADC_fit_pretty", 900, 800 );
    can_ADC_fit_pretty->SetBottomMargin(0.15);
    can_ADC_fit_pretty->SetLeftMargin(0.15);
    can_ADC_fit_pretty->SetTopMargin(0.05);
    can_ADC_fit_pretty->SetRightMargin(0.05);

    vec_ADC_pads_merge.resize(N_TRD);
    vec_ADC_pads_merge_corr.resize(N_TRD);
    tg_mean_vs_pressure.resize(N_TRD);
    h_slopes.resize(N_TRD);
    h_slopes_over_b.resize(N_TRD);

    vec_h_fit_params.resize(2);
    vec_h_fit_params[0] = new TProfile( "vec_h_fit_params_mean", "vec_h_fit_params_mean", N_TRD, 0, N_TRD);
    vec_h_fit_params[1] = new TProfile( "vec_h_fit_params_sigma", "vec_h_fit_params_sigma", N_TRD, 0, N_TRD);

    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        Int_t sector = i_det/30;

        if(sector != sector_use && sector_use > -1 ) continue;
        cout << i_det << endl;

        vec_ADC_pads_merge[i_det].resize(N_rows_merge);
        vec_ADC_pads_merge_corr[i_det].resize(N_rows_merge);
        tg_mean_vs_pressure[i_det].resize(N_rows_merge);
        //printf("sector: %d, sector_use: %d, N_rows_merge: %d \n",sector,sector_use,N_rows_merge);

        HistName = "h_slopes_over_b_";
        HistName += i_det;
        h_slopes_over_b[i_det] = new TH1F(HistName.Data(),HistName.Data(),70,-0.0012,-0.0004);

        HistName = "h_slopes_";
        HistName += i_det;
        h_slopes[i_det] = new TH1F(HistName.Data(),HistName.Data(),70,-20.0,-5.0);

        printf("Open histograms for detector: %d \n",i_det);
        for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
        {
            vec_ADC_pads_merge[i_det][i_row_merge].resize(N_columns_merge);
            vec_ADC_pads_merge_corr[i_det][i_row_merge].resize(N_columns_merge);
            tg_mean_vs_pressure[i_det][i_row_merge].resize(N_columns_merge);
            for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
            {
                // vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge].resize(N_run_ids);
                // vec_ADC_pads_merge_corr[i_det][i_row_merge][i_col_merge].resize(N_run_ids);
                tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge] = new TGraphErrors();

                // for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
                {
                    HistName = "merge_ADC/vec_ADC_pads_merge_det_";
                    HistName += i_det;
                    HistName += "_row_";
                    HistName += i_row_merge;
                    HistName += "_col_";
                    HistName += i_col_merge;
                    // if ( i_det== 1 && i_row_merge==1 && i_col_merge==1 ) 
                    // {
                    //     cout << "Open histogram: " << HistName.Data() << endl;
                    // }
                    vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge] = (TH1F*)inputfile->Get(HistName.Data());
                    
                    // if ( i_det== 0 && i_row_merge==0 && i_col_merge==0 ) 
                    // {
                    //     if( !vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge] ) continue;
                    //     cout << "Open histogram: " << HistName.Data() << endl;
                    //     vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge]->Draw();
                    //     return;
                    // }
                }
            }
        }
    }

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void init_ADC_pad_wise(Int_t sector_use, Int_t mode_use, TString outfilename)
{
    cout << "" << endl;
    printf("init_ADC_pad_wise for sector: %d \n",sector_use);


    // TString file_Krypton_name = main_data_dir;
    //TString file_Krypton_name = "Krypton_sector_";
    //file_Krypton_name += sector_use;
    //file_Krypton_name += "_mode_";
    //file_Krypton_name += mode_use;
    //file_Krypton_name += ".root";

    TString file_Krypton_name = main_data_dir;
    file_Krypton_name += outfilename;
    file_Krypton = new TFile(file_Krypton_name.Data(),"RECREATE");



    nt_Krypton = new TNtuple("nt_Krypton","nt_Krypton","det:col:row:chi2:mean:sigma:amplitude");

    h2D_ADC_mean_main_Krypton_xy_TRD          = new TH2F("h2D_ADC_mean_main_Krypton_xy_TRD","h2D_ADC_mean_main_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    h2D_ADC_sigma_main_Krypton_xy_TRD         = new TH2F("h2D_ADC_sigma_main_Krypton_xy_TRD","h2D_ADC_sigma_main_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    h2D_ADC_amplitude_main_Krypton_xy_TRD     = new TH2F("h2D_ADC_amplitude_main_Krypton_xy_TRD","h2D_ADC_amplitude_main_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    h2D_ADC_chi2_main_Krypton_xy_TRD          = new TH2F("h2D_ADC_chi2_main_Krypton_xy_TRD","h2D_ADC_chi2_main_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);

    // h2D_ADC_mean_second_Krypton_xy_TRD        = new TH2F("h2D_ADC_mean_second_Krypton_xy_TRD","h2D_ADC_mean_second_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    // h2D_ADC_sigma_second_Krypton_xy_TRD       = new TH2F("h2D_ADC_sigma_second_Krypton_xy_TRD","h2D_ADC_sigma_second_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    // h2D_ADC_amplitude_second_Krypton_xy_TRD   = new TH2F("h2D_ADC_amplitude_second_Krypton_xy_TRD","h2D_ADC_amplitude_second_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    // h2D_ADC_chi2_second_Krypton_xy_TRD        = new TH2F("h2D_ADC_chi2_second_Krypton_xy_TRD","h2D_ADC_chi2_second_Krypton_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);

    HistName = "can_ADC_sector_";
    HistName += sector_use;
    can_ADC = new TCanvas(HistName.Data(),HistName.Data(),10,10,680,600);




    vec_ADC_pads.resize(N_TRD);
    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        // printf("Open single pad ADC spectra, i_det: %d \n",i_det);
        Int_t sector = i_det/30;
        if(sector != sector_use && sector_use > -1 ) continue;
        // cout << i_det << endl;

        //printf("Create single pad ADC histogram, i_det: %d \n",i_det);
        vec_ADC_pads[i_det].resize(N_TRD_rows);
        for(Int_t i_row = 0; i_row < N_TRD_rows; i_row++)
        {
            vec_ADC_pads[i_det][i_row].resize(N_TRD_columns);
            for(Int_t i_col = 0; i_col < N_TRD_columns; i_col++)
            {
                HistName = "pad_ADC/";
                HistName += "vec_ADC_pads_det_";
                HistName += i_det;
                HistName += "_row_";
                HistName += i_row;
                HistName += "_col_";
                HistName += i_col;
                //printf("i_row: %d, i_col: %d \n",i_row,i_col);
                //vec_ADC_pads[i_det][i_row][i_col] = (TH1F*) inputfile->Get( HistName.Data() ); // inefficient
                if(vec_ADC_pads[i_det][i_row][i_col])
                {
                    vec_ADC_pads[i_det][i_row][i_col]->Sumw2();
                    vec_ADC_pads[i_det][i_row][i_col]->SetLineColor(kBlack);
                    vec_ADC_pads[i_det][i_row][i_col]->SetLineWidth(2);
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(0.0,6500.0);
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetTitleOffset(1.0);
                    vec_ADC_pads[i_det][i_row][i_col]->GetYaxis()->SetTitleOffset(1.6);
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetLabelSize(0.06);
                    vec_ADC_pads[i_det][i_row][i_col]->GetYaxis()->SetLabelSize(0.06);
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetTitleSize(0.06);
                    vec_ADC_pads[i_det][i_row][i_col]->GetYaxis()->SetTitleSize(0.06);
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetNdivisions(505,'N');
                    vec_ADC_pads[i_det][i_row][i_col]->GetYaxis()->SetNdivisions(505,'N');
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->CenterTitle();
                    vec_ADC_pads[i_det][i_row][i_col]->GetYaxis()->CenterTitle();
                    vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetTitle("ADC");
                    vec_ADC_pads[i_det][i_row][i_col]->GetYaxis()->SetTitle("counts");
                    // vec_ADC_pads[i_det][i_row][i_col]->Rebin(2);
                    // if( i_det == 0 ) {
                    //     printf( "\nhistogram: %s\n", HistName.Data() );
                    // }
                }
            }
        }

    }



    //--------------------------------------
    // Efficient read in of histograms
    inputfile->cd("pad_ADC/");
    inputfile->ReadAll("dirs*"); //read all objects in current directory in memory
    TList* list = inputfile->GetList();

    TIter next(inputfile->GetList());
    TObject *obj;
    while(( obj = next()))
    {
        if(obj->IsFolder())
        {
            printf("list \n");
            TDirectoryFile* dir = (TDirectoryFile*)obj;
            cout << "Nkeys: " << dir ->GetNkeys() << endl;
            Int_t Nkeys = dir ->GetNkeys();
            dir ->ReadAll("dirs*");
            TList* list_dir = dir->GetListOfKeys();
            //list_dir->Print();

            TIter next_list(dir->GetListOfKeys());
            TObject *obj_list;
            Int_t i_key = 0;
            while(( obj_list = next_list()))
            {
                // if(i_key % 1000 == 0) cout << "." << flush;
                // if(i_key % 10000 == 0) printf("i_key %d out of %d \n",i_key,Nkeys);
                //cout << obj_list ->GetName() << ", " << obj_list->ClassName() << endl;
                TString obj_name = (TString)obj_list ->GetName();
                obj_name.Replace(0,17,"");
                //cout << obj_name.Data() << ", first _: " << obj_name.First("_") << endl;
                TString sub_det = obj_name(0,obj_name.First("_"));
                //cout << sub_det.Data() << endl;
                Int_t string_det = sub_det.Atoi();

                obj_name.Replace(0,obj_name.First("_")+1,"");
                obj_name.Replace(0,obj_name.First("_")+1,"");
                //cout << obj_name.Data() << ", first _: " << obj_name.First("_") << endl;
                TString sub_row = obj_name(0,obj_name.First("_"));
                //cout << sub_row.Data() << endl;
                Int_t string_row = sub_row.Atoi();

                obj_name.Replace(0,obj_name.First("_")+1,"");
                obj_name.Replace(0,obj_name.First("_")+1,"");
                //cout << obj_name.Data() << ", first _: " << obj_name.First("_") << endl;
                Int_t string_col = obj_name.Atoi();


                TKey* key = (TKey*)obj_list;
                //TH1 *h1 = (TH1*)key->ReadObj();
                vec_ADC_pads[string_det][string_row][string_col] = (TH1F*)((TH1F*)key->ReadObj())->Clone();

                Int_t entries = vec_ADC_pads[string_det][string_row][string_col] ->GetEntries();
                //printf("det: %d, row: %d, col: %d, entries: %d \n",string_det,string_row,string_col,entries);
                i_key++;
            }
        }
    }
    //--------------------------------------


}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Merge_ADC_pad_wise(Int_t sector_use, Int_t mode_use)
{
    cout << "Merge_ADC_pad_wise" << endl;

    Int_t counter_hist = 0;
    h_ADC_full_merge_det.resize(N_TRD);
    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {

        //-----------------------
        // Check if this chamber was flagged as bad, e.g. reduced or off HV
        Int_t flag_bad = 0;
        for(Int_t i_bad = 0; i_bad < (Int_t)vec_bad_chambers.size(); i_bad++)
        {
            Int_t i_det_bad = vec_bad_chambers[i_bad];
            if(i_det == i_det_bad)
            {
                flag_bad = 1;
                break;
            }
        }
        if(flag_bad) continue;
        //-----------------------


        //-----------------------
        // Check if this chamber was flagged as bad, e.g. reduced or off HV
        Int_t flag_add_bad = 0;
        for(Int_t i_bad = 0; i_bad < (Int_t)vec_add_bad_chambers.size(); i_bad++)
        {
            Int_t i_det_bad = vec_add_bad_chambers[i_bad];
            if(i_det == i_det_bad)
            {
                flag_add_bad = 1;
                break;
            }
        }
        if(flag_add_bad) continue;
        //-----------------------


        Int_t sector = i_det/30;
        if(sector != sector_use) continue;

        printf("Merge single pad ADC spectra, i_det: %d \n",i_det);
        for(Int_t i_row = 0; i_row < N_TRD_rows; i_row++)
        {
            for(Int_t i_col = 0; i_col < N_TRD_columns; i_col++)
            {
                if(vec_ADC_pads[i_det][i_row][i_col])
                {
                    if(!h_ADC_full_merge_det[i_det])
                    {
                        HistName = "h_ADC_full_merge_det_";
                        HistName += i_det;
                        h_ADC_full_merge_det[i_det] = (TH1F*)vec_ADC_pads[i_det][i_row][i_col]->Clone(HistName.Data());
                    }
                    else
                    {
                        h_ADC_full_merge_det[i_det]->Add(vec_ADC_pads[i_det][i_row][i_col]);
                    }

                    if(counter_hist == 0)
                    {
                        h_ADC_full_merge = (TH1F*)vec_ADC_pads[i_det][i_row][i_col]->Clone("h_ADC_full_merge");
                    }
                    else
                    {
                        h_ADC_full_merge->Add(vec_ADC_pads[i_det][i_row][i_col]);
                    }
                    counter_hist++;
                }
            }
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void init_peak_positions(TString peak_pos_file)
{
    cout << "init_peak_positions" << endl;
    // HistName = main_data_dir;
    //HistName = "fits.root";
    HistName = peak_pos_file;
    TFile* inputfile_peak_positions = TFile::Open(HistName.Data());

    vec_h_fit_params_merged.resize(2);
    vec_h_fit_params_merged[0] = (TProfile*)inputfile_peak_positions->Get("vec_h_fit_params_mean");
    vec_h_fit_params_merged[1] = (TProfile*)inputfile_peak_positions->Get("vec_h_fit_params_sigma");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fit_ADC_pad_wise(Int_t sector_use)
{

    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        Int_t sector = i_det/30;

        if(sector != sector_use && sector_use > -1 ) continue;

        Int_t stack  = (i_det%30)/6;
        Int_t layer  = i_det%6;

        double mean_main  = fabs( vec_h_fit_params_merged[0]->GetBinContent( i_det ) );
        double sigma_main = fabs( vec_h_fit_params_merged[1]->GetBinContent( i_det ) );

        if(sigma_main > 500.0) sigma_main = 200.0;
        if(sigma_main < 50.0)  sigma_main = 200.0;

        double mean_global  = vec_h_fit_params_merged[0]->GetMean(2);
        double sigma_global = vec_h_fit_params_merged[1]->GetMean(2);

        if(mean_global <= 0) mean_global = mean_main;
        if(sigma_global > 500.0 || sigma_global <= 50.0) sigma_global = sigma_main;

        printf("Fit single pad ADC spectra, i_det: %d, mean: %4.3f, sigma: %4.3f \n",i_det,mean_main,sigma_main);


        for(Int_t i_row = 0; i_row < N_TRD_rows; i_row++)
        {
            for(Int_t i_col = 0; i_col < N_TRD_columns; i_col++)
            {
                if( !vec_ADC_pads[i_det][i_row][i_col] ) continue;
                Int_t entries = vec_ADC_pads[i_det][i_row][i_col]->GetEntries();
                if(entries < 20) continue;
                if(entries < 50)  vec_ADC_pads[i_det][i_row][i_col]->Rebin(5);
                // if(entries >= 50) vec_ADC_pads[i_det][i_row][i_col]->Rebin(4);

                // Int_t det = (layer + stack * 6 + sector * 6 * 5);
                Int_t x_TRD         = i_col + sector*144; // 144*18
                Int_t y_TRD         = i_row + stack*16 + layer*16*5;


                //-------------------------------
                // Fit main Krypton peak

                Int_t    max_bin  = vec_ADC_pads[i_det][i_row][i_col]->GetMaximumBin();
                Double_t bin_cont = vec_ADC_pads[i_det][i_row][i_col]->GetBinContent(max_bin);
                Double_t bin_cent = vec_ADC_pads[i_det][i_row][i_col]->GetBinCenter(max_bin);

                Double_t pos_left_FWHM = 0.0;
                for(Int_t i_bin = max_bin; i_bin >= 1; i_bin--)
                {
                    Double_t bin_cont_left = vec_ADC_pads[i_det][i_row][i_col]->GetBinContent(i_bin);
                    if(bin_cont_left <= 0.5*bin_cont)
                    {
                        pos_left_FWHM = vec_ADC_pads[i_det][i_row][i_col]->GetBinCenter(i_bin);
                        break;
                    }
                }
                Double_t pos_right_FWHM = 0.0;
                for(Int_t i_bin = max_bin; i_bin < vec_ADC_pads[i_det][i_row][i_col]->GetNbinsX(); i_bin++)
                {
                    Double_t bin_cont_right = vec_ADC_pads[i_det][i_row][i_col]->GetBinContent(i_bin);
                    if(bin_cont_right <= 0.5*bin_cont)
                    {
                        pos_right_FWHM = vec_ADC_pads[i_det][i_row][i_col]->GetBinCenter(i_bin);
                        break;
                    }
                }

                Double_t FWHM = (pos_right_FWHM - pos_left_FWHM)*0.5;
                Double_t fit_start = bin_cent-2.0*FWHM;
                Double_t fit_stop  = bin_cent+2.5*FWHM;
                Int_t flag_low_ADC_peak = 0;
                if(bin_cent < 1200.0 && FWHM > 200.0)
                {
                    FWHM = 150.0;
                    fit_start = bin_cent-1.0*FWHM;
                    flag_low_ADC_peak = 1;
                }

                //vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(mean_main-1.5*sigma_main,mean_main+2.0*sigma_main);
                //vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(bin_cent-1.5*sigma_main,bin_cent+2.0*sigma_main);
                vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(bin_cent-2.5*FWHM,bin_cent+2.5*FWHM);
                Double_t rms = vec_ADC_pads[i_det][i_row][i_col]->GetRMS();


                for(Int_t i = 0; i < 4; i++)
                {
                    func_Gauss_fit->ReleaseParameter(i);
                    func_Gauss_fit->SetParameter(i,0.0);
                    func_Gauss_fit->SetParError(i,0.0);
                }
                func_Gauss_fit->SetParameter(0,bin_cont);
                func_Gauss_fit->SetParLimits(0,0,10000);
                func_Gauss_fit->SetParameter(1,bin_cent);
                //func_Gauss_fit->SetParameter(2,sigma_main);
                //func_Gauss_fit->SetParameter(2,rms*0.8);
                func_Gauss_fit->SetParameter(2,FWHM);
                func_Gauss_fit->FixParameter(3,0.0);

                //printf("i_det, i_row, i_col: {%d,%d,%d}, set params-> amplitude, mean, sigma: {%f,%f,%f} \n",i_det,i_row,i_col,bin_cont,bin_cent,sigma_main);

                vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(0.0,15000.0);
                //vec_ADC_pads[i_det][i_row][i_col]->Fit("func_Gauss_fit","QMN","",bin_cent-1.0*sigma_main,bin_cent+1.5*sigma_main);
                vec_ADC_pads[i_det][i_row][i_col]->Fit("func_Gauss_fit","QMN","",fit_start,fit_stop);

                double amplitude = func_Gauss_fit->GetParameter(0);
                double mean      = func_Gauss_fit->GetParameter(1);
                double sigma     = fabs(func_Gauss_fit->GetParameter(2));

                if(i_det == 5) printf("   x/y TRD: {%d, %d}, row, col: {%d, %d}, entries: %d, rms: %4.3f, FWHM: %4.3f, bin_cent(max): %4.3f, bin_cont: %4.3f, amp, mean, sigma: {%4.3f, %4.3f, %4.3f}, range: {%4.3f, %4.3f} \n",x_TRD,y_TRD,i_row,i_col,entries,rms,FWHM,bin_cent,bin_cont,amplitude,mean,sigma,fit_start,fit_stop);

                // //printf("i_det, i_row, i_col: {%d,%d,%d}, amplitude, mean, sigma: {%f,%f,%f} \n",i_det,i_row,i_col,amplitude,mean,sigma);

                // // Limit fit parameters to something useful, this is only needed if the fit fails
                if ( mean > 8000. || mean < 200. ) mean = mean_global; //2000.;
                if ( sigma > 1000.0 || sigma < 50.0 ) sigma = sigma_global; //300;
                if ( amplitude <= 0. ) amplitude = bin_cont;

                // for(Int_t i = 0; i < 4; i++)
                // {
                //     func_Gauss_fit->SetParameter(i,0.0);
                //     func_Gauss_fit->SetParError(i,0.0);
                // }

                func_Gauss_fit->SetParameter(0,amplitude);
                func_Gauss_fit->SetParLimits(0,0,10000);
                func_Gauss_fit->SetParameter(1,mean);
                func_Gauss_fit->SetParameter(2,sigma);

                for(Int_t i = 0; i < 3; i++)
                {
                    func_Gauss_fit->SetParError(i,0.0);
                }

                Double_t left_sigma = 2.5;
                if(flag_low_ADC_peak)
                {
                    left_sigma = 1.0;
                }

                func_Gauss_fit->SetRange(mean-left_sigma*sigma,mean+2.5*sigma);
                vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(mean-left_sigma*sigma,mean+2.5*sigma);
                vec_ADC_pads[i_det][i_row][i_col]->Fit("func_Gauss_fit","QMN","",mean-left_sigma*sigma,mean+2.5*sigma);

                amplitude = func_Gauss_fit->GetParameter(0);
                mean      = func_Gauss_fit->GetParameter(1);
                sigma     = fabs(func_Gauss_fit->GetParameter(2));

                double chi2      = func_Gauss_fit->GetChisquare();
                int    NDF       = func_Gauss_fit->GetNDF();

                vec_ADC_pads[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(0.0,15000.0);
                //if(i_det == 0 && i_row == 0 && i_col == 10) printf("i_det, i_row, i_col: {%d,%d,%d}, entries: %d, amplitude, mean, sigma, chi2, NDF: {%f, %f, %f, %f, %d} \n",i_det,i_row,i_col,entries,amplitude,mean,sigma,chi2,NDF);

                if(NDF > 0) chi2 /= (Double_t)NDF;
                else chi2 = -1.0;

                if(i_det == 5) printf("  --> row, col: {%d, %d}, entries: %d, mean: %4.3f, sigma: %4.3f \n",i_row,i_col,entries,mean,sigma);

                if(mean > 200.0 && mean < 8000.0 && sigma > 20.0 && sigma < 1000.0 && amplitude > 0.0)
                {
                    h2D_ADC_mean_main_Krypton_xy_TRD     ->SetBinContent(x_TRD+1,y_TRD+1,mean);
                    h2D_ADC_sigma_main_Krypton_xy_TRD    ->SetBinContent(x_TRD+1,y_TRD+1,sigma);
                    h2D_ADC_amplitude_main_Krypton_xy_TRD->SetBinContent(x_TRD+1,y_TRD+1,amplitude);
                    h2D_ADC_chi2_main_Krypton_xy_TRD     ->SetBinContent(x_TRD+1,y_TRD+1,chi2);
                }
                else
                {
                    printf("\033[0;31m");
                    printf("   WARNING: not filled! det: %d, entries: %d \n",i_det,entries);
                    printf("   x/y TRD: {%d, %d}, row, col: {%d, %d}, entries: %d, rms: %4.3f, bin_cent(max): %4.3f, bin_cont: %4.3f, amp, mean, sigma: {%4.3f, %4.3f, %4.3f} \n",x_TRD,y_TRD,i_row,i_col,entries,rms,bin_cent,bin_cont,amplitude,mean,sigma);

                    printf("\033[0m");
                }

                if(i_det == 23 && i_row == 13 && i_col == 104)
                {
                    can_ADC->cd()->SetLeftMargin(0.2);
                    can_ADC->cd()->SetRightMargin(0.15);
                    can_ADC->cd()->SetBottomMargin(0.2);
                    can_ADC->cd();
                    vec_ADC_pads[i_det][i_row][i_col]->DrawCopy("h");
                    func_Gauss_fit->SetLineColor(kRed);
                    func_Gauss_fit->SetLineStyle(1);
                    func_Gauss_fit->SetLineWidth(2);
                    func_Gauss_fit->DrawClone("same");
                }
                if(i_det == 23 && i_row == 13 && i_col == 107)
                {
                    can_ADC->cd()->SetLeftMargin(0.2);
                    can_ADC->cd()->SetRightMargin(0.15);
                    can_ADC->cd()->SetBottomMargin(0.2);
                    can_ADC->cd();
                    vec_ADC_pads[i_det][i_row][i_col]->SetLineColor(kBlue);
                    vec_ADC_pads[i_det][i_row][i_col]->DrawCopy("same h");
                    func_Gauss_fit->SetLineColor(kBlue);
                    func_Gauss_fit->SetLineStyle(1);
                    func_Gauss_fit->SetLineWidth(2);
                    func_Gauss_fit->DrawClone("same");
                }
                // //-------------------------------

                // //printf("i_det: %d \n",i_det);
                nt_Krypton->Fill((Float_t)i_det, (Float_t)i_col, (Float_t)i_row, (Float_t)chi2, (Float_t)mean, (Float_t)sigma, (Float_t)amplitude);

            } // i_col
        } // i_row

        // if ( i_det == 10 ) break;
    } // i_det
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void save_Krypton_calibration(Int_t sector_use, Int_t mode_use)
{
    cout << "save_Krypton_calibration" << endl;

    file_Krypton->cd();
    cout << "meow" << endl;

    if(mode_use == 1) {


        h2D_ADC_mean_main_Krypton_xy_TRD->Write();
        h2D_ADC_sigma_main_Krypton_xy_TRD->Write();
        h2D_ADC_amplitude_main_Krypton_xy_TRD->Write();
        h2D_ADC_chi2_main_Krypton_xy_TRD->Write();

        can_ADC->Write();

        nt_Krypton->Write();


        file_Krypton->cd();
        file_Krypton->mkdir("histos");
        file_Krypton->cd("histos");

        for(Int_t i_col = 0; i_col < 144; i_col++) {
            for(int i_row = 0; i_row < 16; i_row++ ) {
                if(vec_ADC_pads[52][i_row][i_col]) vec_ADC_pads[52][i_row][i_col]->Write();
            }
        }

        // for(Int_t i_col = 0; i_col < 144; i_col++)
        // {
        //     if(vec_ADC_pads[0][0][i_col])
        //     {
        //         vec_ADC_pads[0][0][i_col] ->Write();
        //     }

        //     if(vec_ADC_pads[5][5][i_col])
        //     {
        //         vec_ADC_pads[5][5][i_col] ->Write();
        //     }
        // }
        file_Krypton->cd();
    }

    // if(mode_use == 2)
    // {
    //     h_ADC_full_merge->Write();

    //     for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    //     {
    //         Int_t sector = i_det/30;
    //         if(sector != sector_use) continue;
    //         if(h_ADC_full_merge_det[i_det]) h_ADC_full_merge_det[i_det]->Write();
    //     }
    // }

    file_Krypton->Close();
    printf("file_Krypton closed \n");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fit_ADC_pad_merge(Int_t sector_use)
{
    cout << "" << endl;
    cout << "fit_ADC_pad_merge started" << endl;
    Double_t ADC_start_value = 500.0; // search range for main Krypton peak
    Double_t ADC_stop_value  = 8000.0;
    Int_t    ADC_bin_start = 1;
    Int_t    ADC_bin_stop  = -1;
    for(Int_t i_det = 0; i_det < N_TRD; i_det++) {

        Int_t sector = i_det/30;
        Int_t stack = (i_det%30)/6;
        Int_t layer = i_det%6;

        if(sector != sector_use && sector_use > -1 ) continue;
        cout << "i_det: " << i_det << endl;

        vec_merge_fit_par[i_det].resize(N_rows_merge);
        
        for(Int_t i_row = 0; i_row < N_rows_merge; i_row++) {

            //printf("i_row: %d \n",i_row);
            vec_merge_fit_par[i_det][i_row].resize(N_columns_merge);
            for(Int_t i_col = 0; i_col < N_columns_merge; i_col++) {

                //printf("i_col: %d \n",i_col);
                vec_merge_fit_par[i_det][i_row][i_col].resize(3); // Gauss fit parameters
                if(!vec_ADC_pads_merge[i_det][i_row][i_col]) continue;
                // cout << "meow" << endl;
                Double_t N_entries = vec_ADC_pads_merge[i_det][i_row][i_col]->GetEntries();
                if(vec_ADC_pads_merge[i_det][i_row][i_col]->GetEntries() < 100) continue;
                vec_ADC_pads_merge[i_det][i_row][i_col]->Rebin(8);
                vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(ADC_start_value,ADC_stop_value);
                ADC_bin_start = vec_ADC_pads_merge[i_det][i_row][i_col]->FindBin(ADC_start_value);
                ADC_bin_stop  = vec_ADC_pads_merge[i_det][i_row][i_col]->FindBin(ADC_stop_value);
                Double_t integral = vec_ADC_pads_merge[i_det][i_row][i_col]->Integral(ADC_bin_start,ADC_bin_stop);
                if(integral < 100) continue;
                Int_t bin_max    = vec_ADC_pads_merge[i_det][i_row][i_col]->GetMaximumBin();
                Double_t ADC_max = vec_ADC_pads_merge[i_det][i_row][i_col]->GetBinCenter(bin_max);
                Double_t Amp_max = vec_ADC_pads_merge[i_det][i_row][i_col]->GetBinContent(bin_max);
                double sigma_start = 600;

                //--------------------------------------
                for(Int_t i = 0; i < 4; i++)
                {
                    func_Gauss_fit->SetParameter(i,0.0);
                    func_Gauss_fit->SetParError(i,0.0);
                }
                func_Gauss_fit->SetParameter(0,Amp_max);
                func_Gauss_fit->SetParLimits(0, 0, Amp_max*3);
                func_Gauss_fit->SetParameter(1,ADC_max);
                // func_Gauss_fit->SetParLimits(1, 0, 1e4);
                func_Gauss_fit->SetParameter(2,sigma_start);
                // func_Gauss_fit->SetParLimits(2,10, 1e3);
                func_Gauss_fit->FixParameter(3,0.0);

                // if ( i_det == 301 ) {
                //     vec_ADC_pads_merge[i_det][i_row][i_col]->Fit("func_Gauss_fit","MN","",ADC_max-600.0,ADC_max+800.0);
                // } else {
                //     vec_ADC_pads_merge[i_det][i_row][i_col]->Fit("func_Gauss_fit","QMN","",ADC_max-600.0,ADC_max+800.0);
                // }
                vec_ADC_pads_merge[i_det][i_row][i_col]->Fit("func_Gauss_fit","QMN","",ADC_max-600.0,ADC_max+600.0);

            
                Double_t fit_par[3] = {0.0};
                for(Int_t i = 0; i < 3; i++)
                {
                    fit_par[i] = func_Gauss_fit->GetParameter(i);
                }

                if ( i_det == 301 ) printf("det, row, col: {%d, %d, %d}, N_entries: %4.3f, bin_max: %d, ADC_max: %4.3f, Amp_max: %4.3f, amp: %4.3f, mean: %4.3f, sigma: %4.3f,  \n",i_det,i_row,i_col,N_entries,bin_max,ADC_max,Amp_max,fit_par[0],fit_par[1],fit_par[2]);


                // // Limit fit parameters to something useful, this is only needed if the fit fails
                if ( fit_par[1] > ADC_stop_value || fit_par[1] < ADC_start_value ) fit_par[1] = 2000.;
                if ( fit_par[2] > 800.0 || fit_par[2] < 50.0 ) fit_par[2] = 300;

                for(Int_t i = 0; i < 3; i++) {
                    func_Gauss_fit->SetParameter(i,fit_par[i]);
                    func_Gauss_fit->SetParError(i,0.0);
                }
                func_Gauss_fit->SetParLimits(0,0,Amp_max*3);
                // func_Gauss_fit->SetParLimits(1, 0, 1e4);
                // func_Gauss_fit->SetParLimits(2, 100, 1e3);
                // func_Gauss_fit->SetParLimits(1, fit_par[1]-5*fit_par[2], fit_par[1]+5*fit_par[2]);

                vec_ADC_pads_merge[i_det][i_row][i_col]->Fit("func_Gauss_fit","QMN","",fit_par[1]-1*fit_par[2],fit_par[1]+1.*fit_par[2]);

                Double_t fit_par_err[3] = {0.0};
                for(Int_t i = 0; i < 3; i++)
                {
                    fit_par[i]     = func_Gauss_fit->GetParameter(i);
                    fit_par_err[i] = func_Gauss_fit->GetParError(i);
                    vec_merge_fit_par[i_det][i_row][i_col][i] = fit_par[i];
                }


                printf("   ----> mean value: %4.3f \n",fit_par[1]);

                if ( i_det == 0 ) {
                    can_ADC_fit->cd(i_col+N_columns_merge*i_row + 1);
                    can_ADC_fit->cd(i_col+N_columns_merge*i_row + 1)->SetRightMargin(0.1);
                    can_ADC_fit->cd(i_col+N_columns_merge*i_row + 1)->SetLeftMargin(0.1);
                    can_ADC_fit->cd(i_col+N_columns_merge*i_row + 1)->SetTopMargin(0.05);
                    can_ADC_fit->cd(i_col+N_columns_merge*i_row + 1)->SetBottomMargin(0.15);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->SetRangeUser(0.0,3950.0);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->SetStats(0);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->SetTitle("");
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->SetTitleOffset(1.2);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetYaxis()->SetTitleOffset(1.2);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->SetLabelSize(0.06);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetYaxis()->SetLabelSize(0.06);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->SetTitleSize(0.06);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetYaxis()->SetTitleSize(0.06);
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->SetNdivisions(505,'N');
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetYaxis()->SetNdivisions(505,'N');
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetXaxis()->CenterTitle();
                    vec_ADC_pads_merge[i_det][i_row][i_col]->GetYaxis()->CenterTitle();
                    vec_ADC_pads_merge[i_det][i_row][i_col]->DrawCopy("hist");
                    func_Gauss_fit->SetLineColor(kRed);
                    func_Gauss_fit->SetLineStyle(1);
                    func_Gauss_fit->SetLineWidth(3);
                    func_Gauss_fit->DrawCopy("same");

                }

                // if ( i_det == 0 && i_row == 0 && i_col == 0 ) {
                //     can_ADC_fit_pretty->cd();
                //     vec_ADC_pads_merge[i_det][i_row][i_col]->DrawCopy("hist");
                //     func_Gauss_fit->DrawCopy("same");

                //     TPaveText *pt = new TPaveText(.65,.5,.95,.9);
                //     // pt->SetFillStyle(0);
                //     pt->AddText( Form("det: %i", i_det) );
                //     pt->AddText( Form("row: %i", i_row) );
                //     pt->AddText( Form("column: %i", i_col) );
                //     pt->AddText( "" );
                //     pt->AddText( Form("mean: %.1f", func_Gauss_fit->GetParameter(1) ) );
                //     pt->AddText( Form("sigma: %.1f", func_Gauss_fit->GetParameter(2) ) );
                //     pt->Draw();

                // }

                Int_t x_bin = i_col + N_columns_merge*sector + 1;
                Int_t y_bin = i_row + N_rows_merge*layer + N_TRD_layers*N_rows_merge*stack + 1;
                if( fit_par[1] >= 0.0 && fit_par[1] < 10000.0 && fit_par[2] >= 0.0 && fit_par[2] < 1000 )
                {
                    // cout << "setting h2D_ADC_vs_merged_pads" << endl;
                    // Int_t x_bin = i_col + N_columns_merge*sector + 1;
                    // Int_t y_bin = i_row + N_rows_merge*layer + N_TRD_layers*N_rows_merge*stack + 1;
                    h2D_ADC_vs_merged_pads->SetBinContent(x_bin,y_bin,fit_par[1]);
                    // cout << "set" << endl;
                    vec_h_fit_params[0]->Fill( i_det, fit_par[1] );
                    vec_h_fit_params[1]->Fill( i_det, fit_par[2] );
                }
                h2D_ADC_vs_merged_pads_all->SetBinContent(x_bin,y_bin,fit_par[1]);
                //--------------------------------------

                double amplitude = func_Gauss_fit->GetParameter(0);
                double mean      = func_Gauss_fit->GetParameter(1);
                double sigma     = fabs(func_Gauss_fit->GetParameter(2));

                double chi2      = func_Gauss_fit->GetChisquare();
                int    NDF       = func_Gauss_fit->GetNDF();

                if(NDF > 0) {
                    chi2 /= (double)NDF;
                } else {
                    chi2 = -1.0;
                }

                if ( chi2 > 12. ) {
                    chi2 = -1.0;
                    // amplitude = -999.;
                }

                nt_Krypton->Fill((Float_t)i_det, (Float_t)i_col, (Float_t)i_row, (Float_t)chi2, (Float_t)mean, (Float_t)sigma, (Float_t)amplitude);


            } // col

            // printf("\n\n Done fitting.\n");

        } // row
    
    } // det
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fit_mean_vs_pressure(Int_t sector_use)
{
    cout << "fit_mean_vs_pressure started" << endl;
    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        Int_t sector = i_det/30;
        if(sector != sector_use) continue;
        vec_merge_mean_pressure_fit_par[i_det].resize(N_rows_merge);
        for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
        {
            vec_merge_mean_pressure_fit_par[i_det][i_row_merge].resize(N_columns_merge);
            for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
            {
                vec_merge_mean_pressure_fit_par[i_det][i_row_merge][i_col_merge].resize(2);
                for(Int_t i = 0; i < 6; i++)
                {
                    func_Poly_fit->SetParameter(i,0.0);
                    func_Poly_fit->SetParError(i,0.0);
                    if(i > 1) func_Poly_fit->FixParameter(i,0.0);
                }
                func_Poly_fit->SetParameter(0,15000.0);
                func_Poly_fit->SetParameter(1,-12.0);

                tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge]->Fit("func_Poly_fit","QMN","",920.0,1020.0);

                for(Int_t i = 0; i < 2; i++)
                {
                    vec_merge_mean_pressure_fit_par[i_det][i_row_merge][i_col_merge][i] = func_Poly_fit->GetParameter(i);
                }

                Int_t bin_number = i_col_merge + N_columns_merge*i_row_merge + N_rows_merge*N_columns_merge*i_det + 1;
                h_ADC_vs_pressure_fit_par0->SetBinContent(bin_number,vec_merge_mean_pressure_fit_par[i_det][i_row_merge][i_col_merge][0]);
                h_ADC_vs_pressure_fit_par1->SetBinContent(bin_number,vec_merge_mean_pressure_fit_par[i_det][i_row_merge][i_col_merge][1]);
            }
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void plot_ADC_pad_merge(Int_t det_full)
{
    cout << "plot_ADC_pad_merge started" << endl;
    vector<TCanvas*> vec_can_TRD_pads_merge;
    vec_can_TRD_pads_merge.resize(N_run_ids);
    for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
    {
        HistName = "vec_can_TRD_pads_merge_det";
        HistName += i_run_ids;
        vec_can_TRD_pads_merge[i_run_ids] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1820,1120);
        vec_can_TRD_pads_merge[i_run_ids]->SetRightMargin(0.05);
        vec_can_TRD_pads_merge[i_run_ids]->SetTopMargin(0.02);
        vec_can_TRD_pads_merge[i_run_ids]->SetLeftMargin(0.2);
        vec_can_TRD_pads_merge[i_run_ids]->SetBottomMargin(0.18);
        vec_can_TRD_pads_merge[i_run_ids]->Divide(N_columns_merge,N_rows_merge,0.0,0.0);

        if(i_run_ids == 0)
        {
            vec_ADC_pads_merge_sum_run_ids.resize(N_rows_merge);
            vec_ADC_pads_merge_sum_run_ids_corr.resize(N_rows_merge);
        }

        for(Int_t i_row = 0; i_row < N_rows_merge; i_row++)
        {
            if(i_run_ids == 0)
            {
                vec_ADC_pads_merge_sum_run_ids[i_row].resize(N_columns_merge);
                vec_ADC_pads_merge_sum_run_ids_corr[i_row].resize(N_columns_merge);
            }
            for(Int_t i_col = 0; i_col < N_columns_merge; i_col++)
            {
                Int_t iPad = i_col + N_columns_merge*i_row + 1;
                vec_can_TRD_pads_merge[i_run_ids]->cd(iPad);
                //vec_can_TRD_pads_merge[i_run_ids]->cd(iPad)->SetRightMargin(0.0);
                //vec_can_TRD_pads_merge[i_run_ids]->cd(iPad)->SetTopMargin(0.0);
                //vec_can_TRD_pads_merge[i_run_ids]->cd(iPad)->SetLeftMargin(0.0);
                //vec_can_TRD_pads_merge[i_run_ids]->cd(iPad)->SetBottomMargin(0.0);

                if(!vec_ADC_pads_merge[det_full][i_row][i_col]) continue;
                Int_t max_bin = vec_ADC_pads_merge[det_full][i_row][i_col]->GetMaximumBin();
                Double_t bin_max = vec_ADC_pads_merge[det_full][i_row][i_col]->GetBinContent(max_bin);

                HistName = "vec_ADC_pads_merge_sum_run_ids_r";
                HistName += i_row;
                HistName += "_c";
                HistName += i_col;
                if(i_run_ids == 0)
                {
                    vec_ADC_pads_merge_sum_run_ids[i_row][i_col]      = (TH1F*)vec_ADC_pads_merge[det_full][i_row][i_col]->Clone(HistName.Data());
                }
                else vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->Add(vec_ADC_pads_merge[det_full][i_row][i_col],1.0);

                vec_ADC_pads_merge[det_full][i_row][i_col]->SetLineColor(kBlack);
                vec_ADC_pads_merge[det_full][i_row][i_col]->SetLineWidth(2);
                vec_ADC_pads_merge[det_full][i_row][i_col]->SetFillStyle(3001);
                vec_ADC_pads_merge[det_full][i_row][i_col]->SetFillColor(kGray);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->SetRangeUser(0.0,6500.0);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->SetRangeUser(0.0,bin_max*1.2);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->SetTitleOffset(1.0);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->SetTitleOffset(1.0);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->SetLabelSize(0.06);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->SetLabelSize(0.09);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->SetTitleSize(0.06);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->SetTitleSize(0.06);
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->SetNdivisions(505,'N');
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->SetNdivisions(505,'N');
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->CenterTitle();
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->CenterTitle();
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetXaxis()->SetTitle("");
                vec_ADC_pads_merge[det_full][i_row][i_col]->GetYaxis()->SetTitle("");
                vec_ADC_pads_merge[det_full][i_row][i_col]->DrawCopy("h");


                HistName = vec_ADC_pads_merge[det_full][i_row][i_col]->GetName();
                HistName += "_corr";
                vec_ADC_pads_merge_corr[det_full][i_row][i_col][i_run_ids] = (TH1F*)vec_ADC_pads_merge[det_full][i_row][i_col]->Clone(HistName.Data());
                vec_ADC_pads_merge_corr[det_full][i_row][i_col][i_run_ids]->Reset();

                HistName = "Det: ";
                HistName += det_full;
                HistName += ", row: ";
                HistName += i_row;
                HistName += ", col: ";
                HistName += i_col;
                HistName += ", run idx: ";
                HistName += i_run_ids;
                plotTopLegend((char*)HistName.Data(),0.95,0.91,0.081,kBlack,0.0,42,1,31); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

                Double_t pressure = vec_pressure[i_run_ids];
                HistName = "";
                sprintf(NoP,"%2.1f",pressure);
                HistName += NoP;
                HistName += " hPa";
                plotTopLegend((char*)HistName.Data(),0.95,0.82,0.081,kBlack,0.0,42,1,31); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


                //--------------------------------------
                for(Int_t i = 0; i < 3; i++)
                {
                    func_Gauss_fit->SetParameter(i,vec_merge_fit_par[det_full][i_row][i_col][i]);
                }

                Double_t mean  = vec_merge_fit_par[det_full][i_row][i_col][1];
                Double_t sigma = vec_merge_fit_par[det_full][i_row][i_col][2];
                func_Gauss_fit->SetRange(mean-1.5*sigma,mean+2.5*sigma);

                func_Gauss_fit->SetLineColor(kRed);
                func_Gauss_fit->SetLineStyle(1);
                func_Gauss_fit->SetLineWidth(2);
                func_Gauss_fit->DrawCopy("same");
                //--------------------------------------

                HistName = "#mu = ";
                sprintf(NoP,"%2.1f",mean);
                HistName += NoP;
                plotTopLegend((char*)HistName.Data(),0.95,0.73,0.081,kBlack,0.0,42,1,31); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1



            }
        }

        vec_can_TRD_pads_merge[i_run_ids]->cd(68);
        plotTopLegend((char*)"ADC",0.42,0.03,0.17,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        vec_can_TRD_pads_merge[i_run_ids]->cd(28);
        plotTopLegend((char*)"counts",0.1,0.3,0.2,kBlack,90.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }
    cout << "Plot plot_ADC_pad_merge finished" << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void plot_mean_vs_pressure(Int_t sector_use)
{
    vector<TLegend*> leg_mean_vs_pressure;
    vec_c_mean_vs_pressure.resize(N_TRD_layers*N_TRD_stacks);
    leg_mean_vs_pressure.resize(N_TRD_layers*N_TRD_stacks);

    Int_t det_counter = 0;
    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        Int_t sector = i_det/30;
        if(sector != sector_use) continue;

        HistName = "vec_c_mean_vs_pressure_";
        HistName += det_counter;
        vec_c_mean_vs_pressure[det_counter] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1820,1120);
        vec_c_mean_vs_pressure[det_counter]->Divide(3,3,0.01,0.01);

        // Get max/min values in both directions
        Double_t min_max_x_y[2][2] = {0.0};
        Int_t counter = 0;
        for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
        {
            for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
            {
                Double_t x_val, y_val;
                for(Int_t i_point = 0; i_point < tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge]->GetN(); i_point++)
                {
                    tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge]->GetPoint(i_point,x_val,y_val);
                    if(counter == 0)
                    {
                        min_max_x_y[0][0] = x_val;
                        min_max_x_y[1][0] = x_val;
                        min_max_x_y[0][1] = y_val;
                        min_max_x_y[1][1] = y_val;
                    }
                    else
                    {
                        if(x_val < min_max_x_y[0][0]) min_max_x_y[0][0] = x_val;
                        if(x_val > min_max_x_y[1][0]) min_max_x_y[1][0] = x_val;
                        if(y_val < min_max_x_y[0][1]) min_max_x_y[0][1] = y_val;
                        if(y_val > min_max_x_y[1][1]) min_max_x_y[1][1] = y_val;
                    }
                    counter++;
                }
            }
        }


        TH1F* h_dummy_mean_vs_pressure = new TH1F("h_dummy_mean_vs_pressure","h_dummy_mean_vs_pressure",400,920,1020);
        h_dummy_mean_vs_pressure->GetXaxis()->SetRangeUser(0.0,6500.0);
        h_dummy_mean_vs_pressure->GetXaxis()->SetTitleOffset(1.0);
        h_dummy_mean_vs_pressure->GetYaxis()->SetTitleOffset(1.0);
        h_dummy_mean_vs_pressure->GetXaxis()->SetLabelSize(0.06);
        h_dummy_mean_vs_pressure->GetYaxis()->SetLabelSize(0.06);
        h_dummy_mean_vs_pressure->GetXaxis()->SetTitleSize(0.06);
        h_dummy_mean_vs_pressure->GetYaxis()->SetTitleSize(0.06);
        h_dummy_mean_vs_pressure->GetXaxis()->SetNdivisions(505,'N');
        h_dummy_mean_vs_pressure->GetYaxis()->SetNdivisions(505,'N');
        h_dummy_mean_vs_pressure->GetXaxis()->CenterTitle();
        h_dummy_mean_vs_pressure->GetYaxis()->CenterTitle();
        h_dummy_mean_vs_pressure->GetXaxis()->SetTitle("pressure [hPa]");
        h_dummy_mean_vs_pressure->GetYaxis()->SetTitle("<ADC>");
        h_dummy_mean_vs_pressure->GetXaxis()->SetRangeUser(min_max_x_y[0][0]-5.0,min_max_x_y[1][0]+5.0);
        h_dummy_mean_vs_pressure->GetYaxis()->SetRangeUser(min_max_x_y[0][1]*0.95,min_max_x_y[1][1]*1.05);

        leg_mean_vs_pressure[det_counter] = new TLegend(0.69,0.62,0.88,0.90); // x1,y1,x2,y2
        leg_mean_vs_pressure[det_counter]->SetBorderSize(0);
        leg_mean_vs_pressure[det_counter]->SetFillColor(0);
        leg_mean_vs_pressure[det_counter]->SetTextSize(0.045);

        for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
        {
            vec_c_mean_vs_pressure[det_counter]->cd(i_col_merge+1)->SetRightMargin(0.085);
            vec_c_mean_vs_pressure[det_counter]->cd(i_col_merge+1)->SetTopMargin(0.065);
            vec_c_mean_vs_pressure[det_counter]->cd(i_col_merge+1)->SetLeftMargin(0.11);
            vec_c_mean_vs_pressure[det_counter]->cd(i_col_merge+1)->SetBottomMargin(0.12);

            h_dummy_mean_vs_pressure->DrawCopy("h");

            for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
            {
                tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge]->SetMarkerStyle(20);
                tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge]->SetMarkerColor(arr_color_row_merge[i_row_merge]);
                tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge]->DrawClone("same p");

                if(i_col_merge == 0)
                {
                    HistName = "row merge: ";
                    HistName += i_row_merge;
                    leg_mean_vs_pressure[det_counter]->AddEntry(tg_mean_vs_pressure[i_det][i_row_merge][i_col_merge],HistName.Data(),"p");
                }

                for(Int_t i = 0; i < 6; i++)
                {
                    func_Poly_fit->SetParameter(i,0.0);
                    func_Poly_fit->SetParError(i,0.0);
                    if(i < 2)
                    {
                        func_Poly_fit->SetParameter(i,vec_merge_mean_pressure_fit_par[i_det][i_row_merge][i_col_merge][i]);
                    }
                }

                func_Poly_fit->SetLineColor(arr_color_row_merge[i_row_merge]);
                func_Poly_fit->SetLineStyle(1);
                func_Poly_fit->SetLineWidth(2);
                func_Poly_fit->DrawCopy("same");
            }

            HistName = "Det: ";
            HistName += i_det;
            HistName += ", col merge: ";
            HistName += i_col_merge;
            plotTopLegend((char*)HistName.Data(),0.25,0.95,0.062,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }

        vec_c_mean_vs_pressure[det_counter]->cd(1);
        leg_mean_vs_pressure[det_counter]->Draw();
        h_dummy_mean_vs_pressure->Delete();
        det_counter++;
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void plot_slopes(Int_t sector_use)
{
    vec_c_slopes.resize(2); // slopes, slopes/b

    for(Int_t i_s = 0; i_s < 2; i_s++)
    {
        HistName = "vec_c_slopes_";
        HistName += i_s;
        vec_c_slopes[i_s] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1820,1120);

        vec_c_slopes[i_s]->SetRightMargin(0.1);
        vec_c_slopes[i_s]->SetTopMargin(0.1);
        vec_c_slopes[i_s]->SetLeftMargin(0.15);
        vec_c_slopes[i_s]->SetBottomMargin(0.15);
        vec_c_slopes[i_s]->Divide(N_TRD_layers,N_TRD_stacks);
    }

    Int_t det_counter = 0;
    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        Int_t sector = i_det/30;
        if(sector != sector_use) continue;
        Int_t stack = (i_det%30)/6;
        Int_t layer = i_det%6;

        Int_t iPad = layer + 6*stack + 1;

        //printf("i_det: %d, sector: %d, stack: %d, layer: %d, iPad: %d \n",i_det,sector,stack,layer,iPad);

        for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
        {
            for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
            {

                Double_t fit_par[2] = {0.0};
                for(Int_t i = 0; i < 2; i++)
                {
                    fit_par[i] = vec_merge_mean_pressure_fit_par[i_det][i_row_merge][i_col_merge][i];
                }
                Double_t slope_over_b = 0.0;
                if(fit_par[0] > 0.0) // should be always positive
                {
                    slope_over_b = fit_par[1]/fit_par[0];
                }
                //printf("slope_over_b: %f \n",slope_over_b);
                h_slopes[i_det]     ->Fill(fit_par[1]);
                h_slopes_over_b[i_det]->Fill(slope_over_b);
            }
        }

        //--------------------
        vec_c_slopes[0]->cd(iPad);
        vec_c_slopes[0]->cd(iPad)->SetRightMargin(0.15);
        vec_c_slopes[0]->cd(iPad)->SetTopMargin(0.1);
        vec_c_slopes[0]->cd(iPad)->SetLeftMargin(0.15);
        vec_c_slopes[0]->cd(iPad)->SetBottomMargin(0.15);
        h_slopes[i_det]->SetLineColor(kBlack);
        h_slopes[i_det]->SetLineWidth(2);
        h_slopes[i_det]->SetFillStyle(3001);
        h_slopes[i_det]->SetFillColor(kGray);
        //h_slopes[i_det]->GetXaxis()->SetRangeUser(0.0,6500.0);
        h_slopes[i_det]->GetXaxis()->SetTitleOffset(1.0);
        h_slopes[i_det]->GetYaxis()->SetTitleOffset(1.0);
        h_slopes[i_det]->GetXaxis()->SetLabelSize(0.06);
        h_slopes[i_det]->GetYaxis()->SetLabelSize(0.06);
        h_slopes[i_det]->GetXaxis()->SetTitleSize(0.06);
        h_slopes[i_det]->GetYaxis()->SetTitleSize(0.06);
        h_slopes[i_det]->GetXaxis()->SetNdivisions(505,'N');
        h_slopes[i_det]->GetYaxis()->SetNdivisions(505,'N');
        h_slopes[i_det]->GetXaxis()->CenterTitle();
        h_slopes[i_det]->GetYaxis()->CenterTitle();
        h_slopes[i_det]->GetXaxis()->SetTitle("slope");
        h_slopes[i_det]->GetYaxis()->SetTitle("counts");
        h_slopes[i_det]->DrawCopy("h");

        for(Int_t i = 0; i < 4; i++)
        {
            func_Gauss_fit->SetParameter(i,0.0);
            func_Gauss_fit->SetParError(i,0.0);
        }
        func_Gauss_fit->SetParameter(0,10);
        func_Gauss_fit->SetParameter(1,-13.0);
        func_Gauss_fit->SetParameter(2,3.0);
        func_Gauss_fit->FixParameter(3,0.0);
        h_slopes[i_det]->Fit("func_Gauss_fit","QMN","",-18.0,-9.0);

        Double_t mean  = func_Gauss_fit->GetParameter(1);
        Double_t sigma = fabs(func_Gauss_fit->GetParameter(2));
        func_Gauss_fit->SetRange(mean-1.5*sigma,mean+2.5*sigma);

        func_Gauss_fit->SetLineColor(kRed);
        func_Gauss_fit->SetLineStyle(1);
        func_Gauss_fit->SetLineWidth(2);
        func_Gauss_fit->DrawCopy("same");

        Double_t rel_res = 0.0;
        if(mean != 0.0)
        {
            rel_res = sigma/mean;
        }
        sprintf(NoP,"%2.3f",rel_res);
        HistName = "#sigma/#mu = ";
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.6,0.82,0.05,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        //--------------------



        //--------------------
        vec_c_slopes[1]->cd(iPad);
        vec_c_slopes[1]->cd(iPad)->SetRightMargin(0.05);
        vec_c_slopes[1]->cd(iPad)->SetTopMargin(0.1);
        vec_c_slopes[1]->cd(iPad)->SetLeftMargin(0.13);
        vec_c_slopes[1]->cd(iPad)->SetBottomMargin(0.12);
        h_slopes_over_b[i_det]->SetLineColor(kBlack);
        h_slopes_over_b[i_det]->SetLineWidth(2);
        h_slopes_over_b[i_det]->SetFillStyle(3001);
        h_slopes_over_b[i_det]->SetFillColor(kGray);
        h_slopes_over_b[i_det]->GetXaxis()->SetRangeUser(-1.0,-0.6);
        h_slopes_over_b[i_det]->GetXaxis()->SetTitleOffset(1.0);
        h_slopes_over_b[i_det]->GetYaxis()->SetTitleOffset(1.0);
        h_slopes_over_b[i_det]->GetXaxis()->SetLabelSize(0.06);
        h_slopes_over_b[i_det]->GetYaxis()->SetLabelSize(0.06);
        h_slopes_over_b[i_det]->GetXaxis()->SetTitleSize(0.06);
        h_slopes_over_b[i_det]->GetYaxis()->SetTitleSize(0.06);
        h_slopes_over_b[i_det]->GetXaxis()->SetNdivisions(505,'N');
        h_slopes_over_b[i_det]->GetYaxis()->SetNdivisions(505,'N');
        h_slopes_over_b[i_det]->GetXaxis()->CenterTitle();
        h_slopes_over_b[i_det]->GetYaxis()->CenterTitle();
        h_slopes_over_b[i_det]->GetXaxis()->SetTitle("slope/b");
        h_slopes_over_b[i_det]->GetYaxis()->SetTitle("counts");
        h_slopes_over_b[i_det]->DrawCopy("h");

        for(Int_t i = 0; i < 4; i++)
        {
            func_Gauss_fit->SetParameter(i,0.0);
            func_Gauss_fit->SetParError(i,0.0);
        }
        func_Gauss_fit->SetParameter(0,10);
        func_Gauss_fit->SetParameter(1,-0.0008);
        func_Gauss_fit->SetParameter(2,0.00001);
        func_Gauss_fit->FixParameter(3,0.0);
        h_slopes_over_b[i_det]->Fit("func_Gauss_fit","QMN","",-0.0009,-0.0007);

        mean  = func_Gauss_fit->GetParameter(1);
        sigma = fabs(func_Gauss_fit->GetParameter(2));
        func_Gauss_fit->SetRange(mean-1.5*sigma,mean+2.5*sigma);

        func_Gauss_fit->SetLineColor(kRed);
        func_Gauss_fit->SetLineStyle(1);
        func_Gauss_fit->SetLineWidth(2);
        func_Gauss_fit->DrawCopy("same");

        rel_res = 0.0;
        if(mean != 0.0)
        {
            rel_res = sigma/mean;
        }
        sprintf(NoP,"%2.3f",rel_res);
        HistName = "#sigma/#mu = ";
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.6,0.82,0.05,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        //--------------------


        for(Int_t i_s = 0; i_s < 2; i_s++)
        {
            vec_c_slopes[i_s]->cd(iPad);
            HistName = "Det: ";
            HistName += i_det;
            plotTopLegend((char*)HistName.Data(),0.25,0.94,0.062,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void plot_ADC_TRD_chambers(Int_t sector_plot)
{
    vec_can_TRD_det_ADC_sectors.resize(N_TRD_sectors);
    for(Int_t i_sector = 0; i_sector < N_TRD_sectors; i_sector++)
    {
        if(i_sector != sector_plot && sector_plot != -1) continue;
        HistName = "vec_can_TRD_det_ADC_sectors_sec";
        HistName += i_sector;
        vec_can_TRD_det_ADC_sectors[i_sector] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1820,1120);
        vec_can_TRD_det_ADC_sectors[i_sector]->SetRightMargin(0.15);
        vec_can_TRD_det_ADC_sectors[i_sector]->SetTopMargin(0.03);
        vec_can_TRD_det_ADC_sectors[i_sector]->SetLeftMargin(0.17);
        vec_can_TRD_det_ADC_sectors[i_sector]->SetBottomMargin(0.17);
        vec_can_TRD_det_ADC_sectors[i_sector]->Divide(N_TRD_layers,N_TRD_stacks,0.01,0.01);
        for(Int_t i_layer = 0; i_layer < N_TRD_layers; i_layer++)
        {
            for(Int_t i_stack = 0; i_stack < N_TRD_stacks; i_stack++)
            {
                Int_t i_chamber = i_layer + i_stack*N_TRD_layers;
                vec_can_TRD_det_ADC_sectors[i_sector]->cd(i_chamber+1)->SetRightMargin(0.085);
                vec_can_TRD_det_ADC_sectors[i_sector]->cd(i_chamber+1)->SetTopMargin(0.065);
                vec_can_TRD_det_ADC_sectors[i_sector]->cd(i_chamber+1)->SetLeftMargin(0.11);
                vec_can_TRD_det_ADC_sectors[i_sector]->cd(i_chamber+1)->SetBottomMargin(0.12);

                Int_t i_det = (i_layer + i_stack * 6 + i_sector * 6 * 5);
                vec_h_ADC_TRD_chambers[i_det]->SetLineColor(kBlack);
                vec_h_ADC_TRD_chambers[i_det]->SetLineWidth(2);
                vec_h_ADC_TRD_chambers[i_det]->SetFillStyle(3001);
                vec_h_ADC_TRD_chambers[i_det]->SetFillColor(kGray);
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->SetRangeUser(0.0,6500.0);
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->SetTitleOffset(1.0);
                vec_h_ADC_TRD_chambers[i_det]->GetYaxis()->SetTitleOffset(1.0);
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->SetLabelSize(0.06);
                vec_h_ADC_TRD_chambers[i_det]->GetYaxis()->SetLabelSize(0.06);
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->SetTitleSize(0.06);
                vec_h_ADC_TRD_chambers[i_det]->GetYaxis()->SetTitleSize(0.06);
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->SetNdivisions(505,'N');
                vec_h_ADC_TRD_chambers[i_det]->GetYaxis()->SetNdivisions(505,'N');
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->CenterTitle();
                vec_h_ADC_TRD_chambers[i_det]->GetYaxis()->CenterTitle();
                vec_h_ADC_TRD_chambers[i_det]->GetXaxis()->SetTitle("ADC");
                vec_h_ADC_TRD_chambers[i_det]->GetYaxis()->SetTitle("counts");
                vec_h_ADC_TRD_chambers[i_det]->DrawCopy("h");

                if(i_det == det_full)
                {
#if 0
                    printf("i_chamber: %d \n",i_chamber);
                    for(Int_t i_bin = 1; i_bin < h_corr_ADC_full->GetNbinsX(); i_bin++)
                    {
                        Double_t bin_cont = h_corr_ADC_full->GetBinContent(i_bin);
                        Double_t bin_cent = h_corr_ADC_full->GetBinCenter(i_bin);
                        printf("i_bin: %d, bin_cont: %f, bin_cent: %f \n",i_bin,bin_cont,bin_cent);
                    }
#endif
                    Double_t integral_old  = vec_h_ADC_TRD_chambers[i_det]->Integral(1,-1);
                    Double_t integral_new  = h_corr_ADC_full->Integral(1,-1);
                    Double_t bin_width_old = vec_h_ADC_TRD_chambers[i_det]->GetBinWidth(1);
                    Double_t bin_width_new = h_corr_ADC_full->GetBinWidth(1);

                    //if(integral_new > 0.0) h_corr_ADC_full->Scale(integral_old/integral_new);
                    h_corr_ADC_full->Scale(bin_width_old/bin_width_new);
                    h_corr_ADC_full->SetLineColor(kRed);
                    h_corr_ADC_full->DrawCopy("same h");
                }

                HistName = "Det: ";
                HistName += i_det;
                plotTopLegend((char*)HistName.Data(),0.74,0.85,0.062,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            }
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void plot_ADC_pads_merge_sum_run_ids()
{
    //-------------------------------------------
    cout << "plot_ADC_pads_merge_sum_run_ids started" << endl;
    TCanvas* can_pads_merge_sum_run_ids;

    HistName = "can_pads_merge_sum_run_ids_det";
    can_pads_merge_sum_run_ids = new TCanvas(HistName.Data(),HistName.Data(),10,10,1820,1120);
    can_pads_merge_sum_run_ids->SetRightMargin(0.05);
    can_pads_merge_sum_run_ids->SetTopMargin(0.02);
    can_pads_merge_sum_run_ids->SetLeftMargin(0.2);
    can_pads_merge_sum_run_ids->SetBottomMargin(0.18);
    can_pads_merge_sum_run_ids->Divide(N_columns_merge,N_rows_merge,0.0,0.0);


    for(Int_t i_row = 0; i_row < N_rows_merge; i_row++)
    {
        for(Int_t i_col = 0; i_col < N_columns_merge; i_col++)
        {
            Int_t iPad = i_col + N_columns_merge*i_row + 1;
            can_pads_merge_sum_run_ids->cd(iPad);

            Int_t max_bin    = vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetMaximumBin();
            Double_t bin_max = vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetBinContent(max_bin);

            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->SetLineColor(kBlack);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->SetLineWidth(2);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->SetFillStyle(3001);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->SetFillColor(kGray);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->SetRangeUser(0.0,6500.0);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->SetRangeUser(0.0,bin_max*1.2);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->SetTitleOffset(1.0);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->SetTitleOffset(1.0);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->SetLabelSize(0.06);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->SetLabelSize(0.09);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->SetTitleSize(0.06);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->SetTitleSize(0.06);
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->SetNdivisions(505,'N');
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->SetNdivisions(505,'N');
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->CenterTitle();
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->CenterTitle();
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetXaxis()->SetTitle("");
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->GetYaxis()->SetTitle("");
            vec_ADC_pads_merge_sum_run_ids[i_row][i_col]->DrawCopy("h");

            vec_ADC_pads_merge_sum_run_ids_corr[i_row][i_col]->SetLineColor(kRed);
            vec_ADC_pads_merge_sum_run_ids_corr[i_row][i_col]->SetLineWidth(2);
            vec_ADC_pads_merge_sum_run_ids_corr[i_row][i_col]->DrawCopy("same h");

            HistName = "Det: ";
            HistName += det_full;
            HistName += ", row: ";
            HistName += i_row;
            HistName += ", col: ";
            HistName += i_col;
            plotTopLegend((char*)HistName.Data(),0.95,0.91,0.081,kBlack,0.0,42,1,31); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
    }

    can_pads_merge_sum_run_ids->cd(68);
    plotTopLegend((char*)"ADC",0.42,0.03,0.17,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    can_pads_merge_sum_run_ids->cd(28);
    plotTopLegend((char*)"counts",0.1,0.3,0.2,kBlack,90.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    //-------------------------------------------



    //-------------------------------------------
    TCanvas* can_pads_merge_sum_run_ids_example;
    HistName = "can_pads_merge_sum_run_ids_example_det";
    can_pads_merge_sum_run_ids_example = new TCanvas(HistName.Data(),HistName.Data(),10,10,700,600);
    can_pads_merge_sum_run_ids_example->SetRightMargin(0.12);
    can_pads_merge_sum_run_ids_example->SetTopMargin(0.1);
    can_pads_merge_sum_run_ids_example->SetLeftMargin(0.2);
    can_pads_merge_sum_run_ids_example->SetBottomMargin(0.18);

    Int_t i_row_example = 3;
    Int_t i_col_example = 8;

    vec_ADC_pads_merge_sum_run_ids[i_row_example][i_col_example]->GetXaxis()->CenterTitle();
    vec_ADC_pads_merge_sum_run_ids[i_row_example][i_col_example]->GetYaxis()->CenterTitle();
    vec_ADC_pads_merge_sum_run_ids[i_row_example][i_col_example]->GetXaxis()->SetTitle("ADC");
    vec_ADC_pads_merge_sum_run_ids[i_row_example][i_col_example]->GetYaxis()->SetTitle("counts");

    vec_ADC_pads_merge_sum_run_ids[i_row_example][i_col_example]   ->DrawCopy("h");
    vec_ADC_pads_merge_sum_run_ids_corr[i_row_example][i_col_example]->SetLineWidth(2);
    vec_ADC_pads_merge_sum_run_ids_corr[i_row_example][i_col_example]->SetLineColor(kRed);
    vec_ADC_pads_merge_sum_run_ids_corr[i_row_example][i_col_example]->DrawCopy("same h");

    for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
    {
        vec_ADC_pads_merge[det_full][i_row_example][i_col_example]->DrawCopy("same h");
        vec_ADC_pads_merge_corr[det_full][i_row_example][i_col_example][i_run_ids]->SetLineWidth(2);
        vec_ADC_pads_merge_corr[det_full][i_row_example][i_col_example][i_run_ids]->SetLineColor(kRed);
        vec_ADC_pads_merge_corr[det_full][i_row_example][i_col_example][i_run_ids]->DrawCopy("same h");
    }
    //-------------------------------------------

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void save_pressure_correction(Int_t sector_use, Int_t mode_use, TString soutputfile)
{
    cout << "save_pressure_correction" << endl;
    // TString file_pressure_name = main_data_dir;
    // file_pressure_name += "Pressure_corr_par_sector_";
    // file_pressure_name += sector_use;
    // file_pressure_name += "_mode_";
    // file_pressure_name += mode_use;
    // file_pressure_name += ".root";
    TString file_pressure_name = main_data_dir + soutputfile;
    cout << file_pressure_name << endl;
    // return;
    TFile* file_pressure = new TFile(file_pressure_name.Data(),"RECREATE");
    file_pressure->cd();

// #if 0
//     h_ADC_vs_pressure_fit_par0->Write();
//     h_ADC_vs_pressure_fit_par1->Write();

//     for(Int_t i_s = 0; i_s < 2; i_s++)
//     {
//         vec_c_slopes[i_s]->Write();
//     }
//     for(Int_t i_det = 0; i_det < N_TRD_layers*N_TRD_stacks; i_det++)
//     {
//         vec_c_mean_vs_pressure[i_det]->Write();
//     }
// #endif

    h2D_ADC_vs_merged_pads->Write();
    h2D_ADC_vs_merged_pads_all->Write();
    can_ADC_fit->Write();
    can_ADC_fit_pretty->Write();

    vec_h_fit_params[0]->Write();
    vec_h_fit_params[1]->Write();

    // outputTree->SetDirectory(outputFile)
    nt_Krypton->SetDirectory(file_pressure);
    nt_Krypton->Write();

    file_pressure->cd();
    file_pressure->mkdir("merged_ADC_histos");
    file_pressure->cd("merged_ADC_histos");

    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        Int_t sector = i_det/30;
        if(sector != sector_use && sector_use > -1 ) continue;
        // if(!(i_det == 301)) continue;
        for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
        {
            for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
            {
                for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
                {
                    if(!vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge]) continue;
                    vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge]->Write();
                }
            }
        }
    }

    file_pressure->Close();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Init_input(Int_t sector_use, Int_t mode_use, TString sinputfile)
{
    cout << "Init_input" << endl;

    TString indir = main_data_dir;
    TString full_name;
    /*
    if(mode_use == 0)
    {
        // indir += "Output/";
        // indir += "Output/Merged_ADC_pressure/";
        TString file_name = "KrHistOutput.root";
        // full_name = indir;
        // full_name = file_name;
        full_name = "KrHistOutput_501483_mode0_sector";
        full_name += sector_use;
        full_name += ".root";

        full_name = "/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/output/2021/Output/BKrHistOutput_501483_mode0_sector-1.root";
    }
    if(mode_use == 1)
    {
        // indir += "/Output/Single_Pad/";
        TString file_name = "Merge_Krypton_mode_2.root";
        // full_name = indir;
        // full_name = file_name;
        full_name = "KrHistOutput_501483_mode2_sector";
        full_name += sector_use;
        full_name += ".root";

        full_name = "/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/output/2021/Output/BKrHistOutput_501483_mode2_sector-1.root";
    }
    if(mode_use == 2)
    {
        // indir += "/Output/Single_Pad/";
        //TString file_name = "Merge_Krypton_mode_2.root";
        //TString file_name = "Merge_Krypton_mode_3.root";
        TString file_name = "Merge_Krypton_mode_3_Marco_V2.root";
        // full_name = indir;
        full_name = file_name;
    }
    */
    full_name = main_data_dir;
    full_name += sinputfile;
    printf("\nOpening file: %s\n", full_name.Data() );
    inputfile = TFile::Open(full_name.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Init_pressure_file()
{
    TString pressure_file_name = main_data_dir;
    pressure_file_name += "Output/GRP_output_y2018.root";
    TFile* pressure_file = TFile::Open(pressure_file_name.Data());

    tg_pressure_vs_time       = (TGraph*)pressure_file  ->Get("tg_pressure_vs_time");
    tg_pressure_vs_run_id     = (TGraph*)pressure_file  ->Get("tg_pressure_vs_run_id");

    // Get pressure and plot it
    vec_pressure.resize((Int_t)vec_run_ids.size());
    for(Int_t i_point = 0; i_point < tg_pressure_vs_run_id->GetN(); i_point++)
    {
        Double_t pressure, run_id;
        tg_pressure_vs_run_id->GetPoint(i_point,run_id,pressure);
        //printf("i_point: %d, run_id: %f, pressure: %f \n",i_point,run_id,pressure);
        UInt_t run_id_pressure = (UInt_t)run_id;
        //for(UInt_t run_id_stored : vec_run_ids) // range based loop
        for(Int_t i_run_ids = 0; i_run_ids < (Int_t)vec_run_ids.size(); i_run_ids++)
        {
            UInt_t run_id_stored = vec_run_ids[i_run_ids];
            if(run_id_pressure == run_id_stored)
            {
                vec_pressure[i_run_ids] = pressure;
                break;
            }
        }
    }

    cout << "" << endl;
    cout << "---------------------------------------------" << endl;
    N_run_ids = (Int_t)vec_run_ids.size();
    printf("Total number of run ids: %d \n",N_run_ids);
    for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
    {
        Int_t    run_id_stored = vec_run_ids[i_run_ids];
        Double_t pressure      = vec_pressure[i_run_ids];
        printf("i_run_ids: %d, run_id: %d, pressure: %f \n",i_run_ids,run_id_stored,pressure);
    }
    cout << "---------------------------------------------" << endl;
    cout << "" << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_histo_and_canvas(TH1F* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------