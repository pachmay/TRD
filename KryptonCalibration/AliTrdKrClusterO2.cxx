#include "AliTrdKrClusterO2.h"

ClassImp(AliTrdKrClusterO2);




//------------------------------------------------------------------------------------------------------------------
AliTrdKrClusterO2::AliTrdKrClusterO2()
{

}
//------------------------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------------------------
AliTrdKrClusterO2::~AliTrdKrClusterO2()
{

}
//------------------------------------------------------------------------------------------------------------------


void AliTrdKrClusterO2::createInputList(string sListOfFiles)
{
    ifstream inFileList( sListOfFiles );
    string filename = "";
    while( inFileList >> filename ) {
        cout << filename << endl;
        files_root.push_back( filename );
    }
}

//------------------------------------------------------------------------------------------------------------------
Int_t AliTrdKrClusterO2::createInputList()
{
    // Is creating a list of all root files which are in pinputdir

    printf("AliTrdKrClusterO2::createInputList() \n");

    string dir = string(pinputdir.Data());
    vector<string> files = vector<string>();

    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    if ( sReadFile == "none" ) {
        
        while ((dirp = readdir(dp)) != NULL) {
            files.push_back(string(dirp->d_name));
        }
        closedir(dp);

    } else {
        files.push_back( sReadFile.Data() );
    }

    // while ((dirp = readdir(dp)) != NULL)
    // {
    //     files.push_back(string(dirp->d_name));
    // }
    // closedir(dp);


    Int_t N_files = files.size();
    //printf("N_files: %d \n",N_files);

    string prev_str;
    Int_t good_file_counter = 0;
    std::string str_root = string("root");
    for(Int_t i_file = 0; i_file < N_files; i_file++)
    {
        //cout << files[i_file] << endl;
        Int_t string_length =  files[i_file].length();
        if(string_length > 5)
        {
            std::string str_end = files[i_file].substr(string_length-4,4);
            //cout << str_end << endl;

            if(str_end == str_root)
            {
                files_root.push_back(files[i_file]);
            }
        }
    }

    for(Int_t i_file = 0; i_file < (Int_t)files_root.size(); i_file++)
    {
        cout << "i_file: " << i_file << ", " << files_root[i_file] << endl;
    }

    return 0;
}
//------------------------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------------------------
void AliTrdKrClusterO2::createChain()
{
    printf("AliTrdKrClusterO2::createChain() \n");

    Long64_t entries_save = 0;
    inputChain  = new TChain( sKrClusterTreeName.Data(), sKrClusterTreeName.Data() );
    for(Long64_t i_file = 0; i_file < (Long64_t)files_root.size(); i_file++)
    {
        TString addfile = pinputdir;
        addfile += files_root[i_file];
        //cout << "addfile: " << addfile.Data() << endl;
        inputChain->AddFile(addfile.Data(),-1, sKrClusterTreeName.Data() );
        Long64_t file_entries = inputChain->GetEntries();
        cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
        entries_save = file_entries;
    }
    inputChain->SetBranchAddress( sKrClusterBranchName, &vKrClusterPtr );
    inputChain->SetBranchAddress( sKrClusterTriggerRecordBranchName, &vKrTriggerPtr );

    file_entries_total = inputChain->GetEntries();
    cout << "Total number of events in tree: " << file_entries_total << endl;
 
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void AliTrdKrClusterO2::Init_histograms(Int_t mode)
{
    cout << "Init_histograms" << endl;

    N_columns_merge = N_TRD_columns/merge_N_columns;
    N_rows_merge    = N_TRD_rows/merge_N_rows;

    // printf("N_columns_merge: %d, N_rows_merge: %d \n",N_columns_merge,N_rows_merge);

    h2D_ADC_xy_TRD_merged = new TH2F(
        "h2D_ADC_xy_TRD_merged", "h2D_ADC_xy_TRD_merged",
        N_columns_merge*N_TRD_sectors ,0 ,N_columns_merge*N_TRD_sectors,
        N_TRD_stacks*N_TRD_layers*N_rows_merge, 0, N_TRD_stacks*N_TRD_layers*N_rows_merge
        );
    h2D_ADC_xy_TRD = new TH2F("h2D_ADC_xy_TRD","h2D_ADC_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);

    h2D_xy_TRD_merged= new TH2F( "h2D_xy_TRD_merged", "h2D_xy_TRD_merged",
        N_columns_merge*N_TRD_sectors ,0 ,N_columns_merge*N_TRD_sectors,
        N_TRD_stacks*N_TRD_layers*N_rows_merge, 0, N_TRD_stacks*N_TRD_layers*N_rows_merge
        );

    hProfile2D_ADC_xy_TRD_merged = new TProfile2D(
        "hProfile2D_ADC_xy_TRD_merged", "hProfile2D_ADC_xy_TRD_merged",
        N_columns_merge*N_TRD_sectors ,0 ,N_columns_merge*N_TRD_sectors,
        N_TRD_stacks*N_TRD_layers*N_rows_merge, 0, N_TRD_stacks*N_TRD_layers*N_rows_merge
        );

    // h2D_ADC_pressure = new TH2F( "h2D_ADC_pressure", "h2D_ADC_pressure", 500,0,10000, 160,962,978);

    hSector = new TH1F( "hSector", "hSector", N_TRD_sectors, 0, N_TRD_sectors );

    if(mode == 0 || mode == 1)
    {
        vec_ADC_pads_merge.resize(N_TRD);
        for(Int_t i_det = 0; i_det < N_TRD; i_det++)
        {
            Int_t sector = i_det/30;
            if(sector != useSector && useSector > -1 ) continue;

            //printf("Create merged histograms for pressure correction, det: %d \n",i_det);

            vec_ADC_pads_merge[i_det].resize(N_TRD_rows);
            for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
            {
                vec_ADC_pads_merge[i_det][i_row_merge].resize(N_TRD_columns);
                for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
                {
                    HistName = "vec_ADC_pads_merge_det_";
                    HistName += i_det;
                    HistName += "_row_";
                    HistName += i_row_merge;
                    HistName += "_col_";
                    HistName += i_col_merge;
                    vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge] = new TH1F(HistName.Data(),HistName.Data(),1000,0,15000);
                }
            }
        }
    }

    // if(mode == 1)
    // {
    //     vec_h2D_cls_time_vs_ADC.resize(N_TRD);
    //     vec_h_ADC_cut.resize(N_TRD);
    //     for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    //     {
    //         HistName = "vec_h2D_cls_time_vs_ADC_";
    //         HistName += i_det;
    //         vec_h2D_cls_time_vs_ADC[i_det] = new TH2F(HistName.Data(),HistName.Data(),500,0,15000,30,0,30);
    //     }
    // }


    // if(mode == 3)
    // {
    //     vec_h_ADC_TRD_chambers_runid.resize(N_TRD);
    //     for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    //     {
    //         vec_h_ADC_TRD_chambers_runid[i_det].resize(N_run_ids);
    //         for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
    //         {
    //             HistName = "vec_h_ADC_TRD_chambers_runid_";
    //             HistName += i_det;
    //             HistName += "_runidx_";
    //             HistName += i_run_ids;
    //             vec_h_ADC_TRD_chambers_runid[i_det][i_run_ids] = new TH1F(HistName.Data(),HistName.Data(),1000,0,15000); // 600
    //         }
    //     }


    //     vec_h2D_cls_time_vs_ADC.resize(N_TRD);
    //     vec_h_ADC_cut.resize(N_TRD);
    //     for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    //     {
    //         HistName = "vec_h2D_cls_time_vs_ADC_";
    //         HistName += i_det;
    //         vec_h2D_cls_time_vs_ADC[i_det] = new TH2F(HistName.Data(),HistName.Data(),500,0,15000,30,0,30);


    //         HistName = "vec_h_ADC_cut_";
    //         HistName += i_det;
    //         vec_h_ADC_cut[i_det] = new TH1F(HistName.Data(),HistName.Data(),1000,0,15000);
    //     }


    //     vec_h_ADC_TRD_chambers.resize(N_TRD);
    //     for(Int_t i_hist = 0; i_hist < (Int_t)vec_h_ADC_TRD_chambers.size(); i_hist++)
    //     {
    //         HistName = "vec_h_ADC_TRD_chambers_";
    //         HistName += i_hist;
    //         vec_h_ADC_TRD_chambers[i_hist] = new TH1F(HistName.Data(),HistName.Data(),1000,0,15000); // 600
    //     }

    //     vec_h_ADC_TRD_chambers_cut.resize(N_cuts);
    //     for(Int_t i_hist = 0; i_hist < (Int_t)vec_h_ADC_TRD_chambers_cut.size(); i_hist++)
    //     {
    //         HistName = "vec_h_ADC_TRD_chambers_cut_";
    //         HistName += i_hist;
    //         vec_h_ADC_TRD_chambers_cut[i_hist] = new TH1F(HistName.Data(),HistName.Data(),1000,0,15000); // 600
    //     }

    //     h_det               = new TH1F("h_det","h_det",540,0,540);
    //     // h2D_ADC_xy_TRD      = new TH2F("h2D_ADC_xy_TRD","h2D_ADC_xy_TRD",144*18,0,144*18,5*6*16,0,5*6*16);
    //     h2D_rmsTime_vs_ADC  = new TH2F("h2D_rmsTime_vs_ADC","h2D_rmsTime_vs_ADC",500,0,15000,250,0,5);
    //     h2D_cls_size_vs_ADC = new TH2F("h2D_cls_size_vs_ADC","h2D_cls_size_vs_ADC",500,0,15000,30,0,30);
    //     h2D_cls_dim_vs_ADC  = new TH2F("h2D_cls_dim_vs_ADC","h2D_cls_dim_vs_ADC",500,0,15000,40,0,40);
    //     h2D_cls_time_vs_ADC = new TH2F("h2D_cls_time_vs_ADC","h2D_cls_time_vs_ADC",500,0,15000,30,0,30);
    //     h2D_cls_col_vs_ADC  = new TH2F("h2D_cls_col_vs_ADC","h2D_cls_col_vs_ADC",500,0,15000,15,0,15);
    //     h2D_cls_row_vs_ADC  = new TH2F("h2D_cls_row_vs_ADC","h2D_cls_row_vs_ADC",500,0,15000,15,0,15);
    // }


    if(mode == 2 || mode == 3)
    {
        vec_ADC_pads.resize(N_TRD);
        for(Int_t i_det = 0; i_det < N_TRD; i_det++)
        {
            Int_t sector = i_det/30;
            if(sector != useSector && useSector > -1 ) continue;

            // printf("Create single pad ADC spectra, i_det: %d \n",i_det);
            vec_ADC_pads[i_det].resize(N_TRD_rows);
            for(Int_t i_row = 0; i_row < N_TRD_rows; i_row++)
            {
                // printf("   i_row: %d \n",i_row);
                vec_ADC_pads[i_det][i_row].resize(N_TRD_columns);
                for(Int_t i_col = 0; i_col < N_TRD_columns; i_col++)
                {
                    HistName = "vec_ADC_pads_det_";
                    HistName += i_det;
                    HistName += "_row_";
                    HistName += i_row;
                    HistName += "_col_";
                    HistName += i_col;
                    // printf("      %s\n",HistName.Data());
                    vec_ADC_pads[i_det][i_row][i_col] = new TH1F(HistName.Data(),HistName.Data(),500,0,15000); // 600
                }
            }
        }


        vec_h_deltaTime.resize(2);
        vec_h_deltaColumn.resize(2);
        vec_h_deltaRow.resize(2);
        vec_h_adcRMS.resize(2);
        vec_h_timeRms.resize(2);
        vec_h_adcSum.resize(2);
        vec_h_adcSumOverT.resize(2);

        for(Int_t i_hist = 0; i_hist < 2; i_hist++)
        {
            HistName = "vec_h_deltaTime_";
            HistName += i_hist;
            vec_h_deltaTime[i_hist] = new TH1F(HistName.Data(),HistName.Data(),24,0,24);

            HistName = "vec_h_deltaColumn_";
            HistName += i_hist;
            vec_h_deltaColumn[i_hist] = new TH1F(HistName.Data(),HistName.Data(),10,0,10);

            HistName = "vec_h_deltaRow_";
            HistName += i_hist;
            vec_h_deltaRow[i_hist] = new TH1F(HistName.Data(),HistName.Data(),10,0,10);

            HistName = "vec_h_adcRMS_";
            HistName += i_hist;
            vec_h_adcRMS[i_hist] = new TH1F(HistName.Data(),HistName.Data(),200,0,200);

            HistName = "vec_h_timeRms_";
            HistName += i_hist;
            vec_h_timeRms[i_hist] = new TH1F(HistName.Data(),HistName.Data(),24,0,24);

            HistName = "vec_h_adcSum_";
            HistName += i_hist;
            vec_h_adcSum[i_hist] = new TH1F(HistName.Data(),HistName.Data(),400,0,8000);

            HistName = "vec_h_adcSumOverT_";
            HistName += i_hist;
            vec_h_adcSumOverT[i_hist] = new TH1F(HistName.Data(),HistName.Data(),400,0,8000);
        }

        h_detector_digits = new TH1F("h_detector_digits","h_detector_digits",540,0,540);
    }

    // if(mode == 4)
    // {
    //     vec_ADC_spec_cat.resize(3); // three categories

    //     for(Int_t i_cat = 0; i_cat < 3; i_cat++)
    //     {
    //         vec_ADC_spec_cat[i_cat].resize(N_run_ids);
    //         for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
    //         {
    //             HistName = "vec_ADC_spec_cat_";
    //             HistName += i_cat;
    //             HistName += "_run_id_";
    //             HistName += vec_run_ids[i_run_ids];
    //             vec_ADC_spec_cat[i_cat][i_run_ids] = new TH1F(HistName.Data(),HistName.Data(),500,0,15000);
    //         }
    //     }
    // }
}
//------------------------------------------------------------------------------------------------------------






//------------------------------------------------------------------------------------------------------------------
void AliTrdKrClusterO2::loopClusters(Long64_t nEventsUse, Int_t mode)
{
    if(nEventsUse == -1)                nEventsUse = file_entries_total;
    if(nEventsUse > file_entries_total) nEventsUse = file_entries_total;
    if(nEventsUse < 0)                  nEventsUse = 1;

    vec_h1D_KrClustersDet.resize(N_TRD);
    // vec_h2D_KrClVsTimeDet.resize(N_TRD);
    for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    {
        
        Int_t sector = i_det/30;
        if(sector != useSector && useSector > -1 ) continue;

        HistName = "vec_h1D_KrClustersDet_";
        HistName += i_det;
        vec_h1D_KrClustersDet[i_det] = new TH1F(HistName.Data(),HistName.Data(),500,0,10000);

        // HistName = "vec_h2D_KrClVsTimeDet_";
        // HistName += i_det;
        // vec_h2D_KrClVsTimeDet[i_det] = new TH2F(HistName.Data(),HistName.Data(), 500,0,10000, 160,962,978 );
    }

    Int_t timecut = 28;   // maximum value of delta t
    Int_t columncut = 2;  // maximum value of delta column
    Int_t rowcut = 0;     // maximum value of delta row

    Long64_t nTotalClusters = 0;
    BCData bcData;
    Double_t timeInNs = 0.0;

    printf("\n\n Analyzing run # %i \n Setting initital timestamp: \n", (int)runNumber );
    timeStamp.Print();
    
    Int_t runTimeInSec  = 0;
    Int_t conversion    = 1e9;
    Int_t remainderInNs = 0;
    TTimeStamp clusterStamp;

    for(Long64_t ievent = 0; ievent < nEventsUse; ievent++)
    {
        //cout << "ievent: " << ievent << endl;
        if ( ievent % 1000 == 0 ) printf( "\n Starting event %lld out of %lld \n", ievent,nEventsUse);
        inputChain->GetEntry( ievent );
        Long64_t nNumberKlusters = vKrClusterPtr->size();

     /*   //----------------------------------------
        // Get time information
        clusterStamp = getClusterStamp( vKrTriggerPtr );
        // transform cluster TTimeStamp into TDatime
        int clusterDate = clusterStamp.GetDate();
        int clusterTime = clusterStamp.GetTime(); // gives time down to seconds
        TDatime t( clusterDate, clusterTime );
        double pressure = grP->Eval( t.Convert() ); // get pressure at the cluster time
        // if ( ievent % 1000 == 0 ) printf(" Pressure is %.2f \n", pressure);
    */
        //----------------------------------------
        nTotalClusters += nNumberKlusters;

        for(Long64_t i_cluster = 0; i_cluster < nNumberKlusters; i_cluster++)
        {
            // int sector = detector/30;
            // if(sector != useSector && useSector > -1 ) continue;

            // cout << "i_cluster: " << i_cluster << " out of " << nNumberKlusters << endl;
            Int_t sector = (vKrClusterPtr->at(i_cluster)).getSector();
            if(sector != useSector && useSector > -1 ) continue;

            hSector->Fill( sector );

            Int_t detector    = (vKrClusterPtr->at(i_cluster)).getDetector();
            Int_t column      = (vKrClusterPtr->at(i_cluster)).getColumn();
            Int_t row         = (vKrClusterPtr->at(i_cluster)).getRow();
            Int_t layer       = (vKrClusterPtr->at(i_cluster)).getLayer();
            Int_t stack       = (vKrClusterPtr->at(i_cluster)).getStack();

            Int_t i_row_merge = (Int_t)(row/merge_N_rows);
            Int_t i_col_merge = (Int_t)(column/merge_N_columns);

            Int_t run_id_idx  = 0;

            // if ( ievent % 1000 == 0 ) printf("detector: %d, sector: %d, layer: %d, stack: %d \n",detector,sector,layer,stack);

            Int_t xtrd        = column + sector*144;
            Int_t ytrd        = row + stack*16 + layer*16*5;


            int xtrd_merged = i_col_merge + sector*(N_TRD_columns/merge_N_columns);
            int ytrd_merged = i_row_merge + stack*(N_TRD_rows/merge_N_rows) + layer*(N_TRD_rows/merge_N_rows)*N_TRD_stacks;

            // OK
            Int_t adcRms      = (vKrClusterPtr->at(i_cluster)).getAdcRms();
            Int_t timeRms     = (vKrClusterPtr->at(i_cluster)).getTimeRms();
            Int_t deltaTime   = (vKrClusterPtr->at(i_cluster)).getClSizeTime();
            Int_t deltaRow    = (vKrClusterPtr->at(i_cluster)).getClSizeRow();
            Int_t deltaColumn = (vKrClusterPtr->at(i_cluster)).getClSizeCol();
            Int_t adcMaxA     = (vKrClusterPtr->at(i_cluster)).getAdcMaxA();
            Int_t adcMaxB     = (vKrClusterPtr->at(i_cluster)).getAdcMaxB();
            Int_t adcSum      = (vKrClusterPtr->at(i_cluster)).getAdcSum();
            Int_t adcSumOverT = (vKrClusterPtr->at(i_cluster)).getAdcSumEoverT();

            //if(deltaTime > timecut || deltaColumn > columncut || deltaRow > rowcut)
            if(deltaTime < 28 && deltaRow == 0 && (deltaColumn == 1 || deltaColumn == 2) )
            {
                // if ( ievent % 1000 == 0 ) printf(" filling histos \n");

                //------------------------------------------
                if( mode == 0 )
                {

                    h2D_ADC_xy_TRD_merged->Fill( xtrd_merged, ytrd_merged, adcSumOverT );
                    h2D_xy_TRD_merged->Fill( xtrd_merged, ytrd_merged );
                    h2D_ADC_xy_TRD->Fill( xtrd, ytrd );

                    hProfile2D_ADC_xy_TRD_merged->Fill( xtrd_merged, ytrd_merged, adcSumOverT );

                    vec_h1D_KrClustersDet[detector]->Fill(adcSumOverT);
                    // vec_h2D_KrClVsTimeDet[detector]->Fill(adcSumOverT,pressure);

                    if(i_row_merge >= 0 && i_row_merge < N_rows_merge && i_col_merge >= 0 && i_col_merge < N_columns_merge)
                    {
                        // vec_ADC_pads_merge[detector][i_row_merge][i_col_merge][run_id_idx]->Fill(adcSumOverT);
                        vec_ADC_pads_merge[detector][i_row_merge][i_col_merge]->Fill(adcSumOverT);
                    } else
                    {
                        printf("-----> WARNING: i_row_merge: %d, i_col_merge: %d, row: %d, column: %d \n",i_row_merge,i_col_merge,row,column);
                    }
                } // if mode == 0

                //------------------------------------------
                if(mode == 2)
                {
                    // cout << "detector: " << detector << ", row: " << row << ", column: " << column << ", adcSumOverT: " << adcSumOverT << endl;
                    if(detector >= 0 && detector < N_TRD && row >= 0 && row < N_TRD_rows && column >= 0 && column < N_TRD_columns)
                    {
                        //if(xtrd < 50 && ytrd > 450 && ytrd < 460) printf("detector: %d, adcSumOverT: %d \n",detector,adcSumOverT);
                        if(adcSumOverT > 0) h_detector_digits ->Fill(detector,adcSumOverT);

                        Int_t i_hist = -1;
                        vec_ADC_pads[detector][row][column]->Fill(adcSumOverT);
                        if(detector == 0 && row == 0 && column == 54)
                        {
                            i_hist = 0;
                            //printf("adcSumOverT: %d, adcSum: %d, adcRms: %d, delta(t,r,c): (%d, %d, %d) \n",adcSumOverT,adcSum,adcRms,deltaTime,deltaRow,deltaColumn);
                        }
                        if(detector == 0 && row == 0 && column == 48)
                        {
                            i_hist = 1;
                            //printf("   ----> adcSumOverT: %d, adcSum: %d, adcRms: %d, delta(t,r,c): (%d, %d, %d) \n",adcSumOverT,adcSum,adcRms,deltaTime,deltaRow,deltaColumn);
                        }

                        if(i_hist >= 0)
                        {
                            vec_h_deltaTime[i_hist] ->Fill(deltaTime);
                            vec_h_deltaColumn[i_hist]  ->Fill(deltaColumn);
                            vec_h_deltaRow[i_hist] ->Fill(deltaRow);
                            vec_h_adcRMS[i_hist] ->Fill(adcRms);
                            vec_h_timeRms[i_hist] ->Fill(timeRms);
                            vec_h_adcSum[i_hist] ->Fill(adcSum);
                            vec_h_adcSumOverT[i_hist] ->Fill(adcSumOverT);
                        }
                    }
                } // if mode == 2
            } // if cluster cuts

        } // if nclusters

        // if ( ievent % 1000 == 0 ) printf( "\n Ending event %llu \n", ievent);

    } // if ievent
    printf("nTotalClusters: %lld \n",nTotalClusters);
}
//------------------------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------------------------
void AliTrdKrClusterO2::writeOutput(Int_t mode)
{
    HistName = outputdir;
    if ( sReadFile == "none" )
    {
        HistName += "KrHistOutput_";
        HistName += runNumber;
        HistName += "_mode";
        HistName += mode;
        HistName += "_sector";
        HistName += useSector;
        HistName += ".root";
    }
    else
    {
        HistName += "output_";
        HistName += sReadFile;
    }
    printf("\nWriting into file: %s\n", HistName.Data() );
    TFile* outputfile = new TFile(HistName.Data(),"RECREATE");

    hSector->Write();
    // printf("Write vec_h2D_KrClVsTimeDet \n");
    // for(Int_t i_det = 0; i_det < N_TRD; i_det++)
    // {
    //     vec_h2D_KrClVsTimeDet[i_det]->Write();
    // }

    //------------------------------------------
    cout << "Write histograms" << endl;
    if(mode == 0 || mode == 1) {
        
        outputfile->cd();

        h2D_ADC_xy_TRD_merged->Write();
        h2D_ADC_xy_TRD->Write();
        h2D_xy_TRD_merged->Write();

        hProfile2D_ADC_xy_TRD_merged->Write();
        // h2D_ADC_pressure->Write();

        printf("Write vec_h1D_KrClustersDet \n");
        for(Int_t i_det = 0; i_det < N_TRD; i_det++)
        {
            int sector = i_det/30;
            if( sector != useSector && useSector > -1 ) continue;
            vec_h1D_KrClustersDet[i_det]->Write();
        }

        outputfile->cd();
        outputfile->mkdir("merge_ADC");
        outputfile->cd("merge_ADC");
        printf("Write merged ADC \n");
        
        for(Int_t i_det = 0; i_det < N_TRD; i_det++)
        {
            int sector = i_det/30;
            if( sector != useSector && useSector > -1 ) continue;

            //printf("Save merged histograms for pressure correction, det: %d \n",i_det);
            for(Int_t i_row_merge = 0; i_row_merge < N_rows_merge; i_row_merge++)
            {
                for(Int_t i_col_merge = 0; i_col_merge < N_columns_merge; i_col_merge++)
                {
                    // for(Int_t i_run_ids = 0; i_run_ids < N_run_ids; i_run_ids++)
                    {
                        if(vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge]->GetEntries() <= 0.0) continue;
                        // printf("write, det, row, col, run: {%d,%d,%d,%d} \n",i_det,i_row_merge,i_col_merge,i_run_ids);
                        vec_ADC_pads_merge[i_det][i_row_merge][i_col_merge]->Write();
                    }
                }
            }
        }

        outputfile->cd();
    }
    //------------------------------------------
    if(mode == 2)
    {
        outputfile->cd();

        h_detector_digits ->Write();

        for(Int_t i_hist = 0; i_hist < 2; i_hist++)
        {
            vec_h_deltaTime[i_hist] ->Write();
            vec_h_deltaColumn[i_hist]  ->Write();
            vec_h_deltaRow[i_hist] ->Write();
            vec_h_adcRMS[i_hist] ->Write();
            vec_h_timeRms[i_hist] ->Write();
            vec_h_adcSum[i_hist] ->Write();
            vec_h_adcSumOverT[i_hist] ->Write();
        }


        outputfile->mkdir("pad_ADC");
        outputfile->cd("pad_ADC");

        for(Int_t i_det = 0; i_det < N_TRD; i_det++)
        {
            Int_t sector = i_det/30;
            if(sector != useSector && useSector > -1 ) continue;

            for(Int_t i_row = 0; i_row < N_TRD_rows; i_row++ )
            {
                for(Int_t i_col = 0; i_col < N_TRD_columns; i_col++)
                {
                    if(vec_ADC_pads[i_det][i_row][i_col]->GetEntries() <= 0.0 ) continue;
                    vec_ADC_pads[i_det][i_row][i_col]->Write();
                }
            }
        }
    }

    outputfile->Close();
}
//------------------------------------------------------------------------------------------------------------------




// //------------------------------------------------------------------------------------------------------------------
void AliTrdKrClusterO2::setInitialTimePerRun(int iyear, int imonth, int iday, int ihour, int imin, int isec)
{
    setyear = iyear; 
    setmonth = imonth; 
    setday = iday; 
    sethour = ihour; 
    setmin = imin; 
    setsec = isec; 
    timeStamp.Set( setyear, setmonth, setday, sethour, setmin, setsec, 0, kTRUE, 0);
}
//------------------------------------------------------------------------------------------------------------------




// //------------------------------------------------------------------------------------------------------------------
void AliTrdKrClusterO2::setTreeName(TString treename, TString clustername, TString triggername ) 
{
    sKrClusterTreeName = treename; 
    sKrClusterBranchName = clustername; 
    sKrClusterTriggerRecordBranchName = triggername; 
}
//------------------------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------------------------
TTimeStamp AliTrdKrClusterO2::getClusterStamp(std::vector<o2::trd::KrClusterTriggerRecord>* vTriggerRecordPtr)
{
    int runTimeInSec  = 0;
    int conversion    = 1e9;
    int remainderInNs = 0;
    double timeInNs = 0.0;

    BCData bcData;
    TTimeStamp clusterStamp;

    bcData = vTriggerRecordPtr->at(0).getBCData();  // get bunch crossing data - from there we can get time info
    timeInNs = bcData.bc2ns();  // translate BC data into time in ns
    runTimeInSec = (int)(timeInNs/1000000000);
    remainderInNs = (int)( fmod( timeInNs , conversion ) );
    runTimeInSec += setsec;
    clusterStamp.Set( 2021, 9, setday, sethour, setmin, runTimeInSec, remainderInNs, kTRUE, 0 );

    return clusterStamp;
}
//------------------------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------------------------
// TDatime AliTrdKrClusterO2::convertTimestampToDatime(TTimeStamp timestamp)
// { //// ???
//     TDatime datetime( timestamp.GetDate() , timestamp.GetTime() );
//     return datetime;
// }