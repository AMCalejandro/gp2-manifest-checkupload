try:
    import streamlit as st
    from streamlit.components.v1 import html
    import os
    import sys
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    sys.path.append('utils')
    import generategp2ids
    from customcss import load_css
    from writeread import read_file, to_excel, get_studycode, email_ellie
    from qcutils import detect_multiple_clindups, sample_type_fix
    from plotting import aggridPlotter
except Exception as e:
    print("Some modules are not installed {}".format(e))


os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "secrets/secrets.json"

def jumptwice():
    st.write("##")
    st.write("##")

def app():
    load_css("apps/css/css.css")

    st.markdown("""<div id='link_to_top'></div>""", unsafe_allow_html=True)
    st.markdown('<p class="big-font">GP2 sample manifest self-QC</p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> This is a app tab to self-check the sample manifest and clinical data. </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Download the template from the link below. Once you open the link, go to "File"> "Download" > "xlsx" or "csv" format </p>', unsafe_allow_html=True)
    st.markdown('[Access the data dictionary and templates](https://docs.google.com/spreadsheets/d/1e17_mA8poRSvZCLu9AyBaHIw29uwHHZ0fjgnSddThe0/edit#gid=0)', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Please refer to the second tab (Dictionary) for instructions. </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Please refer to the first tab (sm) to access sample manifest template to fill in. </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Please note all the GP2 required columns must be completed </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Once you have filled in all the columns avavailable in your cohort, please upload the manifest on the side bar to start the QC process </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Once the QC process completes, please move to the Upload tab of the app to store your sample manifest </p>', unsafe_allow_html=True)

    data_file = st.sidebar.file_uploader("Upload Your Sample manifest (CSV/XLSX)", type=['xlsx', 'csv'])
    studycode = get_studycode()
    choice = st.sidebar.selectbox("Genotyping site", 
                                  ["For Fulgent", "For Psomagen", "For NIH", "For LGC", "For UCL", "For DZNE"],
                                  index=None)

    ph_conf = ''
    sex_conf = ''
    race_conf = ''
    fh_conf = ''

    # Columns 
    cols = ['study_type', 'sample_id', 'family_index', 'family_index_relationship', 'sample_type',
            'DNA_volume', 'DNA_conc', 'r260_280',
            'Plate_name', 'Plate_position', 'clinical_id',
            'study_arm', 'diagnosis', 'sex', 'race',
            'age', 'age_of_onset', 'age_at_diagnosis', 'age_at_death', 'age_at_last_follow_up',
            'region', 'comment', 'family_history_pd', 'family_history_other', 'family_history_details', 'alternative_id1', 'alternative_id2']
    required_cols = ['study_type', 'sample_id', 'sample_type', 'clinical_id', 
                     'study_arm', 'diagnosis', 'sex']
    fulgent_cols = ['DNA_volume', 'DNA_conc', 'Plate_name', 'Plate_position']
    #monogenic_cols = ['family_index', 'family_index_relationship', 'age_of_onset', 'family_history_pd', 'family_history_pd_details']
    monogenic_cols = ['age_of_onset', 'family_history_pd']

    # Column values
    gptwo_phenos = ['PD', 'Control', 'Prodromal',
                    'PSP', 'CBD/CBS', 'MSA', 'DLB', 'AD', 'FTD', "VaD", "VaPD"
                    'Population Control', 'Undetermined-MCI', 'Undetermined-Dementia', 'Mix', 'LBD', 'Other']
    allowed_samples = ['Blood (EDTA)', 'Blood (ACD)', 'Blood', 'DNA', 'DNA from Brain',
                        'DNA from blood', 'DNA from FFPE', 'RNA', 'Saliva',
                        'Buccal Swab', 'T-25 Flasks (Amniotic)', 'FFPE Slide',
                        'FFPE Block', 'Fresh tissue', 'Frozen tissue',
                        'Bone Marrow Aspirate', 'Whole BMA', 'CD3+ BMA', 'Other']
    allowed_studyType = ['Case(/Control)', 'Prodromal', 'Genetically Enriched', 'Population Cohort', 'Brain Bank', 'Monogenic']
    allowed_sex = ['Male', 'Female', 'Other/Unknown/Not Reported']
    allowed_race=['American Indian or Alaska Native', 'Asian', 'White', 'Black or African American',
                  'Multi-racial', 'Native Hawaiian or Other Pacific Islander', 'Other', 'Unknown', 'Not Reported']
    allowed_family_history = ['Yes', 'No', 'Not Reported']
    allowed_region_codes=["ABW","AFG","AGO","AIA","ALA","ALB","AND","ARE","ARG","ARM","ASM","ATA","ATF","ATG","AUS","AUT","AZE",
                          "BDI","BEL","BEN","BES","BFA","BGD","BGR","BHR","BHS","BIH","BLM","BLR","BLZ","BMU","BOL","BRA","BRB",
                          "BRN","BTN","BVT","BWA","CAF","CAN","CCK","CHE","CHL","CHN","CIV","CMR","COD","COG","COK","COL","COM",
                          "CPV","CRI","CUB","CUW","CXR","CYM","CYP","CZE","DEU","DJI","DMA","DNK","DOM","DZA","ECU","EGY","ERI",
                          "ESH","ESP","EST","ETH","FIN","FJI","FLK","FRA","FRO","FSM","GAB","GBR","GEO","GGY","GHA","GIB","GIN",
                          "GLP","GMB","GNB","GNQ","GRC","GRD","GRL","GTM","GUF","GUM","GUY","HKG","HMD","HND","HRV","HTI","HUN",
                          "IDN","IMN","IND","IOT","IRL","IRN","IRQ","ISL","ISR","ITA","JAM","JEY","JOR","JPN","KAZ","KEN","KGZ",
                          "KHM","KIR","KNA","KOR","KWT","LAO","LBN","LBR","LBY","LCA","LIE","LKA","LSO","LTU","LUX","LVA","MAC",
                          "MAF","MAR","MCO","MDA","MDG","MDV","MEX","MHL","MKD","MLI","MLT","MMR","MNE","MNG","MNP","MOZ","MRT",
                          "MSR","MTQ","MUS","MWI","MYS","MYT","NAM","NCL","NER","NFK","NGA","NIC","NIU","NLD","NOR","NPL","NRU",
                          "NZL","OMN","PAK","PAN","PCN","PER","PHL","PLW","PNG","POL","PRI","PRK","PRT","PRY","PSE","PYF","QAT",
                          "REU","ROU","RUS","RWA","SAU","SDN","SEN","SGP","SGS","SHN","SJM","SLB","SLE","SLV","SMR","SOM","SPM",
                          "SRB","SSD","STP","SUR","SVK","SVN","SWE","SWZ","SXM","SYC","SYR","TCA","TCD","TGO","THA","TJK","TKL",
                          "TKM","TLS","TON","TTO","TUN","TUR","TUV","TWN","TZA","UGA","UKR","UMI","URY","USA","UZB","VAT","VCT",
                          "VEN","VGB","VIR","VNM","VUT","WLF","WSM","YEM","ZAF","ZMB","ZWE"]

    if data_file is not None:
        
        # Read the data
        jumptwice()
        st.markdown('<p class="big-font">Sample manifest QC </p>', unsafe_allow_html=True)
        df = read_file(data_file)
        df.index = df.index +2
        original_order = df.sample_id.to_list()
        jumptwice()
              
        
        st.subheader("Data columns and sample/clinical IDs check")

        # Select the genotyping site
        jumptwice()
        st.markdown("Adding **Genotyping site** column ")
        if choice==None:
            st.markdown("**Please, choose a Genotyping site on the left bar before starting the columns checks**")
            st.stop()
        else:
            st.markdown("**Genotyping_site** successfully added to the sample manifest")
        
        if choice in ['For Fulgent', 'For Psomagen']:
            required_cols = required_cols + fulgent_cols
        if choice=='For NIH':
            required_cols = required_cols + ['Plate_name', 'Plate_position']
        

        # Check expected columns are present
        missing_cols = np.setdiff1d(cols, df.columns)
        if len(missing_cols)>0:
            st.error(f'{missing_cols} are missing. Please use the template sheet')
            st.stop()
        elif df.shape[1] != len(cols):
            st.error('We have detected more unexpected columns in the input sample manifest')
            st.error(f'{not_required_cols} should not be in the file. Please use the template sheet')
            st.error('Please, try and QC again once this has been sorted')
            not_required_cols = np.setdiff1d(df.columns, cols)
            st.stop()
        else:
            st.markdown('sample manifest **columns** check --> OK')


        # Create genotyping site variable
        df['Genotyping_site'] = choice.replace('For ', '')


        # We do this to detect any weird NA. Then, we efficiently capture them in this block of code and raise error
        df[['sample_id','clinical_id']] = df[['sample_id','clinical_id']].replace('nan', np.nan) # We should detect any weird nas now
        # Check required cols have no na
        df_non_miss_check = df[required_cols].copy()
        if df_non_miss_check.isna().sum().sum()>0:
            st.error('There are some missing entries in the required columns. Please fill the missing cells ')
            st.text('First 30 entries with missing data in any required fields')
            st.write(df_non_miss_check[df_non_miss_check.isna().sum(1)>0].head(30))
            st.stop()
        else:
            st.text('Check no missing data in the required fields --> OK')
        
        # Check monogenic required cols have no na
        if 'Monogenic' in df['study_type'].unique():
            df_mono_non_miss_check = df[df['study_type'] == 'Monogenic'][['sample_id','clinical_id'] + monogenic_cols].copy()
            if df_mono_non_miss_check.isna().sum().sum()>0:
                st.error('There are some missing entries in the required columns for monogenic data. ')
                st.error('Please fill family_history_pd and age_of_onset columns for all monogenic samples in the manifest')
                st.text('First 30 entries with missing data in any required fields')
                st.write(df_mono_non_miss_check[df_mono_non_miss_check.isna().sum(1)>0].head(30))
                st.stop()
            else:
                 st.text('Check no missing data in the required fields for monogenic data --> OK')

        # Transform data type of key columns
        df[['sample_id', 'clinical_id']] = df[['sample_id','clinical_id']].astype(str)
        

        # Check clinical_id data
        jumptwice()
        if df.groupby(['clinical_id']).size().sort_values(ascending=False)[0] > 3:
            detect_multiple_clindups(df)
        dup_ID_this=df.loc[df.duplicated(subset=['clinical_id']), 'clinical_id'].unique()
        if len(dup_ID_this)>0:            
            st.warning(f'Duplicated clinical_id in the manifest: {dup_ID_this}')
            st.warning(f'If this is not expected, please fix it and re upload your sample manifest')
        else:
            st.success('clinical_id in the manifest are all new. No duplication/replication.')


        # Chceck sample_id data
        jumptwice()
        sample_id_dup = df.sample_id[df.sample_id.duplicated()].unique()
        # sample dup check
        if len(sample_id_dup)>0:
            st.text(f'Duplicated sample_id:{sample_id_dup}')
            st.error(f'Unique sample IDs are required (clinical IDs can be duplicated if replicated)')
            st.stop()
        else:
            st.markdown(f'Check there are no **sample_id duplicates** --> OK')
            st.markdown(f'**N** of sample_id (entries):{df.shape[0]}')
            st.markdown(f'**N** of unique clinical_id : {len(df.clinical_id.unique())}')


        # sample type check
        jumptwice()
        st.markdown("**sample_type** check")
        st.write(df.sample_type.astype('str').value_counts())
        not_allowed = np.setdiff1d(df.sample_type.unique(), allowed_samples)
        if len(not_allowed)>0:
            #sample_type_fix(df, allowed_samples)
            sample_type_fix(df = df, allowed_samples = allowed_samples, col = 'sample_type')


        # Study_type check
        jumptwice()
        st.markdown("**study_type** and **study_arm** check")
        st.write(df.groupby(['study_arm', 'study_type']).size().rename('N'))
        st.write(df.pivot_table(index='study_arm', columns='study_type',
                                 values='sample_id', aggfunc='count', margins=True))
        not_allowed = np.setdiff1d(df.study_type.unique(), allowed_studyType)
        if len(not_allowed)>0:
            sample_type_fix(df = df, allowed_samples = allowed_studyType, col = 'study_type')


        # Create study variable
        jumptwice()
        st.markdown("Adding **study** column")
        df['study'] = studycode
        if df['study'].unique() == None:
            st.markdown("**Please, choose a study code on the left bar before assigning GP2 IDs**")
            st.stop()
        else:
            st.markdown("**Study** successfully added to the sample manifest")


        # Adding sampe manifest version using URL data
        if 'manifest_version' in st.experimental_get_query_params():
           mani_vers = 'm' + str(st.experimental_get_query_params()['manifest_version'][0])
           df['manifest_id'] = mani_vers
        else:
            st.error("It seems we cannot assign the version of this manifest. \n Please contact cohort@gp2.org and report this issue")
            st.stop()



        # TEMPORARY CODE
        # DESCRIPTION
        # MONOGENIC SAMPLES ARE ALREADY IN OUR SYSTEM
        # IN ORDER TO INCORPORATE THEM IN THE SYSTEM, WE WILL ACCESS THE JSON, AND WE WILL FILL UP THE SAMPLE MANIFESTS WITH ALREADY ASSIGNED GP2 IDS
        # st.subheader('Find existing GP2 IDs for monogenic samples')
        # if st.session_state['master_get'] == None:
        #     studynames = list(df['study'].unique())
        #     ids_tracker = generategp2ids.master_keyv2(studies = studynames)
        #     study_tracker = ids_tracker[studycode]
        #     st.write(study_tracker)
        #     sampleid_gp2id = {key: value[0] for key, value in study_tracker.items()}
            
        #     df['GP2sampleID'] = df['sample_id'].map(sampleid_gp2id)
        #     df[['GP2ID', 'SampleRepNo']] = df['GP2sampleID'].str.split('_', expand=True)[[0, 2]]
        #     df = df[['GP2sampleID', 'GP2ID', 'SampleRepNo'] + [col for col in df.columns if col not in ['GP2sampleID', 'GP2ID', 'SampleRepNo']]]

        #     missing_gp2ids = df[df['GP2sampleID'].isnull()].copy()
        #     if missing_gp2ids.empty:
        #         st.write("We have been able to find all GP2 IDs for the input sample IDs")
        #         st.write("GP2 IDs assignment --> OK")
        #         st.session_state['df_finalids'] = df
        #         st.session_state['master_get'] = 'DONE'
        #     else:
        #         st.write(" We have detectec a problem looking for GP2 samples IDs for the sample manifest loaded")
        #         st.write("Showing below the samples that should be in our system but are not")
        #         aggridPlotter(missing_gp2ids[['study','study_type','sample_id','clinical_id']])
        #         st.stop()
        
        # else:
        #     df = st.session_state['df_finalids']

        # GP2 IDs assignment
        jumptwice()
        st.subheader('GP2 IDs assignment...')
        if st.session_state['master_get'] == None: # TO ONLY RUN ONCE
            studynames = list(df['study'].unique())
            ids_tracker = generategp2ids.master_keyv2(studies = studynames)
            
            # # Check if this is another QC run of a sample manifest
            all_sids = []
            for v in ids_tracker.values():
                all_sids = all_sids + list(v.keys())
            if len(np.setdiff1d(df['sample_id'].to_list(), all_sids)) == 0:
                st.error("It seems that you are trying to QC the same sample manifest again")
                st.error("Please, contact us on cohort@gp2.org and explain your situation")
                st.stop()
            #     if st.button("Start QC again"):
            #         ids_tracker = generategp2ids.master_remove(studynames, df)
            #         st.experimental_rerun()
            #     st.stop()
            study_subsets = []
            log_new = []
            df['GP2sampleID'] = None
            for study in studynames:
                st.write(f"Getting GP2IDs for {study} samples")
                df_subset = df[df.study==study].copy()
                try:
                    #study_tracker = st.session_state['store_tracker'][study]
                    study_tracker = ids_tracker[study]
                    study_tracker_df = pd.DataFrame.from_dict(study_tracker,
                                                            orient='index',
                                                            columns = ['master_GP2sampleID','clinical_id'])\
                                                    .rename_axis('master_sample_id').reset_index()\
                                                    .astype(str)

                    # Check if any sample ID exists in df_subset.
                    sample_id_unique = pd.merge(study_tracker_df, df_subset,
                                                left_on=['master_sample_id'], right_on=['sample_id'], how='inner')
                    if not sample_id_unique.empty:
                        st.error('We have detected sample ids submitted on previous versions')
                        st.error('Please, correct these sample IDs so that they are unique and resubmit the sample manifest.')
                        st.error('If this is an attempt to re QC a sample manifest, please contact us on cohort@gp2.org')
                        sample_id_unique = sample_id_unique.rename(columns={"clinical_id_y": "clinical_id"})
                        st.dataframe(
                        sample_id_unique[['study','sample_id','clinical_id']].style.set_properties(**{"background-color": "brown", "color": "lawngreen"})
                        )
                        stopapp=True
                    else:
                        stopapp=False
                except:
                    study_tracker = None
                    stopapp = False
                if stopapp:
                    st.stop()

                if bool(study_tracker):
                    # WORK ON DUPLICATED IDS
                    df_subset = df_subset.reset_index()
                    data_duplicated = pd.merge(df_subset, study_tracker_df, on=['clinical_id'], how='inner')
                    df_subset = df_subset.set_index('index')
                    df_subset.index.name = None

                    if data_duplicated.shape[0]>0:
                        new_clinicaldups = True
                        newids_clinicaldups = data_duplicated.groupby('clinical_id')\
                                                        .apply(lambda x: generategp2ids.assign_unique_gp2clinicalids(df_subset,x))

                        if newids_clinicaldups.shape[0]>0:
                            newids_clinicaldups = newids_clinicaldups.reset_index(drop=True)[['study','clinical_id','sample_id','GP2sampleID']]
                            log_new.append(newids_clinicaldups)
                    else:
                        new_clinicaldups = False
                        newids_clinicaldups = pd.DataFrame()

                    # GET GP2 IDs METADATA for new CLINICAL-SAMPLE ID pairs
                    df_newids = df_subset[df_subset['GP2sampleID'].isnull()].reset_index(drop = True).copy()
                    if not df_newids.empty: # Get new GP2 IDs
                        df_wids = df_subset[~df_subset['GP2sampleID'].isnull()].reset_index(drop = True).copy()
                        df_wids['GP2ID'] = df_wids['GP2sampleID'].apply(lambda x: ("_").join(x.split("_")[:-1]))
                        df_wids['SampleRepNo'] = df_wids['GP2sampleID'].apply(lambda x: x.split("_")[-1])#.replace("s",""))

                        n=int(max(study_tracker_df['master_GP2sampleID'].to_list()).split("_")[1])+1
                        df_newids = generategp2ids.getgp2idsv2(df_newids, n, study)
                        df_subset = pd.concat([df_newids, df_wids], axis = 0)
                        study_subsets.append(df_subset)
                        log_new.append(df_newids[['study','clinical_id','sample_id','GP2sampleID']])

                    else: # TO CONSIDER THE CASE IN WHICH WE ONLY HAD DUPLICATE IDS MAPPED ON THE MASTER FILE
                        df_subset['GP2ID'] = df_subset['GP2sampleID'].apply(lambda x: ("_").join(x.split("_")[:-1]))
                        df_subset['SampleRepNo'] = df_subset['GP2sampleID'].apply(lambda x: x.split("_")[-1])#.replace("s",""))
                        study_subsets.append(df_subset)

                # Brand new data - NO STUDY TRACKER FOR THIS COHORT
                else:
                    study = study
                    new_clinicaldups = False # Duplicates from master key json are treated differently to brand new data
                    n = 1
                    df_newids = generategp2ids.getgp2idsv2(df_subset, n, study)
                    study_subsets.append(df_newids)


                # CODE TO UPDATE THE GET FILE WE WILL USE TO UPDATE MASTER JSON
                if (new_clinicaldups) and (newids_clinicaldups.shape[0]>0):
                    tmp = pd.concat([df_newids[['study','clinical_id','sample_id','GP2sampleID']], newids_clinicaldups])
                    tmp['master_value'] = list(zip(tmp['GP2sampleID'],
                                                    tmp['clinical_id']))
                    ids_log = tmp.groupby('study').apply(lambda x: dict(zip(x['sample_id'],
                                                                            x['master_value']))).to_dict()
                else:
                    df_update_master = df_newids.copy()
                    df_update_master['master_value'] = list(zip(df_update_master['GP2sampleID'],
                                                            df_update_master['clinical_id']))
                    ids_log = df_update_master.groupby('study').apply(lambda x: dict(zip(x['sample_id'],
                                                                                    x['master_value']))).to_dict()

                #if st.session_state['master_get'] == None:
                if (isinstance(st.session_state['all_ids'], list)):
                    st.session_state['all_ids'].append( [ids_log, study_tracker] )
                if st.session_state['all_ids'] == None:
                    st.session_state['all_ids'] = [ [ids_log, study_tracker] ]

            # OUT OF FOR LOOP // END OF GP2 IDS ASSIGNMENT. LET'S RESUME df.
            df = pd.concat(study_subsets, axis = 0)
            df = df[list(df)[-3:] + list(df)[:-3]]
            st.write("GP2 IDs assignment... OK")

            #if st.session_state['master_get'] == None:
            if len(log_new) > 0:
                allnew = pd.concat(log_new, axis = 0).reset_index(drop=True)
                st.write("Thanks for uploading a new version of the sample manifest")
                st.write(f'We have detected a total of {allnew.shape[0]} new samples')
                st.write("We have assigned new GP2IDs to those. Showing them below...")
                st.dataframe(
                allnew.style.set_properties(**{"background-color": "brown", "color": "lawngreen"})
                #allnew.style.set_properties(**{"background-color": "brown", "color": "lawngreen"})
                )
            st.session_state['df_finalids'] = df
            st.session_state['master_get'] = 'DONE'

        else:
            df = st.session_state['df_finalids']
            #TODO
            #There is the scenario in which genotyping site or study is changed after assigning GP2 IDs
            # In this scenario, study and genotyping site variables would differ between df and st.session_state['df_finalids']
            # GP2 IDs would refer to the wrong sudy as well
            # We should add some code to consider this scenario


        # Plot the data
        aggridPlotter(df)
        

        # Highlight clinical IDs that were submitted on previous sample manifest
        dup_ID_all=df.loc[df.SampleRepNo!='s1', 'clinical_id'].unique()
        dup_ID_previous = np.setdiff1d(dup_ID_all, dup_ID_this)
        if len(dup_ID_previous)>0:
            st.warning(f'clinical_id previously submitted: {dup_ID_previous}')
            st.warning('If this does not look right, please fix it in your sample manifest, and come back to the app')
        

        # Count by SampleRepNo
        jumptwice()
        st.markdown('Count by **SampleRepNo**: if all clinical IDs are new, all s1')
        st.markdown('Please refer to the dictionary to understand what **SampleRepNo** stands for')
        st.write(df['SampleRepNo'].value_counts())



        ########################
        ### PHENOTYPE FOR QC ###
        ########################
        # Summarise the data
        jumptwice()
        st.subheader('Create Phenotype columns')
        st.text('Show study arm versus diagnosis')
        xtab = df.pivot_table(index='study_arm', columns='diagnosis', margins=True,
                              values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)

        jumptwice()
        st.text('Count per diagnosis')
        st.write(df.diagnosis.astype('str').value_counts())

        # Do the user - gp2 standards mapping
        jumptwice()
        diag = df.diagnosis.dropna().unique()
        n_diag = st.columns(len(diag))
        phenotypes={}
        count_widget = 0
        for i, x in enumerate(n_diag):
            count_widget += 1
            with x:
                mydiag = diag[i]
                diag_index = gptwo_phenos.index(mydiag) if mydiag in gptwo_phenos else None
                phenotypes[mydiag]=x.selectbox(f"[{mydiag}]: For QC, please pick the closest Phenotype",
                                               options = gptwo_phenos,
                                               index = diag_index,
                                               key=count_widget)
                
        df['GP2_phenotype'] = df.diagnosis.map(phenotypes).fillna('Not Assigned') # diagnosis and phenotype relationships are 1:1
        if any(value is not None for value in phenotypes.values()):
            st.text('===  diagnosis x GP2_phenotype ===')
            xtab = df.pivot_table(index='diagnosis', columns='GP2_phenotype', margins=True,
                                  values='sample_id', aggfunc='count', fill_value=0)
            st.write(xtab)
            jumptwice()
            st.text('=== study_arm x study_type x diagnosis x GP2_phenotype===')
            pheno_tab=df.groupby(['study_arm', 'study_type', 'diagnosis', 'GP2_phenotype']).size().rename('N').reset_index()
            st.table(pheno_tab)
        
        # Confirm mapping looks good
        ph_conf = st.checkbox('Confirm Phenotype?')
        if ph_conf:
            if 'Not Assigned' in df.GP2_phenotype.unique():
                st.error('Please assign the phenotype for all the samples')
                st.stop()
            else:
                # Check now if there as some controls with unexpected AAO   
                df_control_aao = df[ (df['GP2_phenotype'].isin(['Control', 'Population Control'])) & (~df['age_of_onset'].isna()) ].copy()
                if not df_control_aao.empty:
                    st.error("We have detected some controls that have age of onset values. This is not possible. Please correct")
                    aggridPlotter(df_control_aao[['study','study_type','sample_id','clinical_id', 'GP2_phenotype', 'age_of_onset']] )
                    st.stop()
                else:
                    st.info('Thank you')

        # Derive phenotype for QC variable
        df['GP2_phenotype_for_qc'] = df['GP2_phenotype'].apply(lambda pheno: pheno if pheno in ['PD','Control'] else 'Other')
        all_studytype = df['study_type'].unique().tolist()
        if 'Genetically Enriched' in all_studytype:
            matching_rows = df[df['study_type'] == 'Genetically Enriched']
            df.loc[matching_rows.index, 'GP2_phenotype_for_qc'] = 'Other'            



        ##################
        ### SEX FOR QC ###
        ##################
        # Summarise the data
        jumptwice()
        st.subheader('Create "biological_sex_for_qc"')
        st.text('Count per sex group')
        st.write(df.sex.astype('str').value_counts())

        # Do the user - GP2 standards mapping
        jumptwice()
        sexes=df.sex.dropna().unique()
        n_sexes = st.columns(len(sexes))
        mapdic={}
        for i, x in enumerate(n_sexes):
            count_widget += 1
            with x:
                sex = sexes[i]
                sex_index = allowed_sex.index(sex) if sex in allowed_sex else None
                mapdic[sex] = x.selectbox(f"[{sex}]: For QC, please pick a word below",
                                        options=allowed_sex, 
                                        index=sex_index, 
                                        key=count_widget)
        
        df['biological_sex_for_qc'] = df['sex'].map(mapdic).fillna('Not Assigned')

        st.text('=== biological_sex_for_qc x sex ===')
        if any(value is not None for value in mapdic.values()):
            xtab = df.pivot_table(index='biological_sex_for_qc', columns='sex', margins=True,
                                    values='sample_id', aggfunc='count', fill_value=0)
            st.write(xtab)

        # Confirm mapping looks good
        sex_conf = st.checkbox('Confirm biological_sex_for_qc?')
        if sex_conf:
            if 'Other/Unknown/Not Reported' in df['biological_sex_for_qc'].unique():
                sex_count = df['biological_sex_for_qc'].value_counts()
                unknown_rate = sex_count['Other/Unknown/Not Reported'] / df.shape[0]
                if unknown_rate > 0.01:
                    st.error("The number of samples with \"Other/Unknown/Not Reported\" sex category is higher than 1% ")
                    st.error("Please check that you selected the right sex values for your samples above ")
                    st.stop()
            elif 'Not Assigned' in df.biological_sex_for_qc.unique():
                st.error('Please assign the sex for all the samples')
                st.stop()
            else:
                st.info('Thank you')



        ###################
        ### RACE FOR QC ###
        ###################
        # Summarise the data
        jumptwice()
        st.subheader('Create "race_for_qc"')
        st.text('Count per race (Not Reported = missing)')
        df['race_for_qc'] = df.race.fillna('Not Reported')
        st.write(df.race.fillna('Not Reported').astype('str').value_counts())

        # Do the user - GP2 standards mapping
        jumptwice()
        races = df.race.dropna().unique()
        nmiss = sum(pd.isna(df.race))
        if nmiss>0:
            st.text(f'{nmiss} entries missing race...')
        mapdic = {'Not Reported':'Not Reported'}
        for race in races:
            count_widget += 1
            race_index=allowed_race.index(race) if race in allowed_race else None
            mapdic[race]=st.selectbox(f"[{race}]: For QC purppose, select the best match from the followings",
                                      options=allowed_race, 
                                      index=race_index, 
                                      key=count_widget)
        
        df['race_for_qc'] = df.race_for_qc.map(mapdic).fillna('Not Assigned')

        st.text('=== race_for_qc X race ===')
        if any(value is not None for value in mapdic.values()):
            dft = df.copy()
            dft['race'] = dft.race.fillna('Not Reported')
            xtab = dft.pivot_table(index='race_for_qc', columns='race', margins=True,
                                values='sample_id', aggfunc='count', fill_value=0)
            st.write(xtab)

        # Confirm mapping looks good
        race_conf = st.checkbox('Confirm race_for_qc?')
        if race_conf:
            if 'Not Assigned' in df.race_for_qc.unique():
                st.error('Please assign the race for all the samples')
                st.stop()
            else:
                st.info('Thank you')



        #################
        ### FH FOR QC ###
        #################
        # Summarise the data
        jumptwice()
        st.subheader('Create "family_history_for_qc"')
        st.text('Count per family_history category (Not Reported = missing)')
        df['family_history_for_qc'] = df.family_history_pd.fillna('Not Reported')
        st.write(df.family_history_pd.fillna('Not Reported').astype('str').value_counts())
        family_historys = df.family_history_pd.dropna().unique()
        nmiss = sum(pd.isna(df.family_history_pd))

        # Do the user - GP2 standards mapping
        jumptwice()
        if nmiss>0:
            st.text(f'{nmiss} entries missing family_history')
        mapdic = {'Not Reported':'Not Reported'}
        
        if len(family_historys)>0:
            n_fhs = st.columns(len(family_historys))
            for i, x in enumerate(n_fhs):
                count_widget += 1
                with x:
                    fh = family_historys[i]
                    fh_index=allowed_family_history.index(fh) if fh in allowed_family_history else None
                    mapdic[fh]=x.selectbox(f'[{fh}]: For QC, family history',
                                           options=allowed_family_history, 
                                           index=fh_index,
                                           key=count_widget)
        
        df['family_history_for_qc'] = df.family_history_for_qc.map(mapdic).fillna('Not Assigned')

        st.text('=== family_history_for_qc X family_history ===')
        if any(value is not None for value in mapdic.values()):
            dft = df.copy()
            dft['family_history_pd'] = dft.family_history_pd.fillna('_Missing')
            xtab = dft.pivot_table(index='family_history_for_qc', columns='family_history_pd', margins=True,
                                    values='sample_id', aggfunc='count', fill_value=0)
            st.write(xtab)

        # Confirm mapping looks good
        fh_conf = st.checkbox('Confirm family_history_for_qc?')
        if fh_conf:
            if 'Not Assigned' in df.family_history_for_qc.unique():
                st.error('Please assign the family for all the samples')
                st.stop()
            else:
                st.info('Thank you')



        #####################
        ### REGION FOR QC ###
        #####################
        # Summarise the data
        jumptwice()
        st.subheader('Create "region_for_qc"')
        st.text('Count per region (Not Reported = missing)')
        df['region_for_qc'] = df.region.fillna('Not Reported')
        st.write(df.region_for_qc.astype('str').value_counts())
        regions = df.region_for_qc.dropna().unique()
        nmiss = sum(pd.isna(df.region))

        # Do the user - GP2 standards mapping
        jumptwice()
        if nmiss>0:
            st.text(f'{nmiss} entries missing for region')
        mapdic = {'Not Reported':'Not Reported'}

        if len(regions)>0:
            st.text('if ISO 3166-3 is available for the region, please provide')
            st.write('https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3')

            if 'Not Reported' in df['region_for_qc'].values:
                st.warning("We have detectec missing values on region column")
                st.warning("For the samples missing the region value, please select the ISO code that corresponds to your STUDY SITE from the URL above")

            n_rgs = st.columns(len(regions))
            for i, x in enumerate(n_rgs):
                count_widget += 1
                with x:
                    region = regions[i]
                    region_index = allowed_region_codes.index(region) if region in allowed_region_codes else None
                    mapdic[region]=x.selectbox(f'[{region}]: For QC, region',
                                           options=allowed_region_codes, 
                                           index=region_index,
                                           key=count_widget)

        df['region_for_qc'] = df.region_for_qc.map(mapdic).fillna('Not Assigned')

        st.text('=== region X  region_for_qc ===')
        if any(value is not None for value in mapdic.values()):
            dft = df.copy()
            dft['region'] = dft.region.fillna('Not Reported')
            xtab = dft.pivot_table(columns='region_for_qc', index='region', margins=True,
                                    values='sample_id', aggfunc='count', fill_value=0)
            st.table(xtab)

        # Confirm mapping looks good
        region_not_well_assigned = np.setdiff1d(df.region_for_qc.unique(), allowed_region_codes) # This should not be necessary anymore, but leave it here for now
        rg_conf = st.checkbox('Confirm region_for_qc?')
        if rg_conf:
            if len(region_not_well_assigned)>0:
                st.error('Please make sure all samples assigned a 3-digit region code. If not available,\
                         please assign the 3-digit code for the principle research site.')
                st.error(f'Need reviews: {region_not_well_assigned}')
                st.stop()
            elif 'Not Assigned' in df.region_for_qc.unique():
                st.error('Please assign the region for all the samples')
                st.stop()
            else:
                st.info('Thank you')



        ##################
        ### PLATE INFO ###
        ##################
        jumptwice()
        st.subheader('Plate Info')
        dft = df.copy()
        dft['Plate_name'] = dft.Plate_name.fillna('_Missing')
        xtab = dft.pivot_table(index='Plate_name',
                            columns='diagnosis', margins=True,
                            values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)

        show_all_duppos = []
        for plate in dft.Plate_name.unique():
            df_plate = dft[dft.Plate_name==plate].copy()
            df_plate_pos = df_plate.Plate_position
            if plate!='_Missing':
                if len(df_plate_pos)>96:
                    st.error('Please make sure, N of samples on plate [{plate}] is =<96')
                    st.stop()
                dup_pos = df_plate_pos[df_plate_pos.duplicated()].unique()
                if len(dup_pos)>0:
                    st.error(f' !!!SERIOUS ERROR!!!  Plate position duplicated position {dup_pos} on plate [{plate}]')
                    show_all_duppos.append((dup_pos, plate))
        if len(show_all_duppos) > 0:
            st.stop()



        ##########################
        ### CHECK NUMERIC COLS ###
        ##########################
        jumptwice()
        st.subheader('Numeric Values')
        numerics_cols = ['DNA_volume', 'DNA_conc', 'r260_280','age', 'age_of_onset', 'age_at_diagnosis', 'age_at_last_follow_up','age_at_death']
        age_cols = ['age', 'age_of_onset', 'age_at_diagnosis']
        for v in numerics_cols:
            if df.dtypes[v] not in ['float64', 'int64']:
                st.error(f'{v} is not numeric')
                st.text("Please, make sure expected numeric columns are stored on a numeric format")
                st.text(f'Expected columns in numeric format - \n * {numerics_cols}')
                st.stop()
            if (df[v] < 0).any():
                st.error(f'{v} has unexpected negative values.  Please correct the values we are showing above')
                aggridPlotter(df[df[v] < 0][ ['study','study_type','sample_id','clinical_id'] + [v] ])
                st.stop()
            if v in age_cols:
                if ((df[v] < 20) | (df[v] > 100)).any():
                    st.error(f'{v} has unexpected high or low values.  Please correct the values we are showing above')
                    aggridPlotter(df[(df[v] < 20) | (df[v] > 100)][ ['study','study_type','sample_id','clinical_id'] + [v] ])
                    st.stop()

        st.text('Numeric chek --> OK.')
        st.text('You can check the distribution with the button below')
        if st.button("Check Distribution"):
            for v in numerics_cols:
                nmiss = df[v].isna().sum()
                vuniq = df[v].dropna().unique()
                nuniq = len(vuniq)
                if nuniq==0:
                    st.text(f'{v} - All missing')
                elif nuniq==1:
                    st.text(f'{v} - One value = {vuniq[0]}, ({nmiss} entries missing)')
                elif nuniq <6:
                    st.write(df[v].value_counts(dropna=False))
                else:
                    st.text(f'{v} - histogram ({nmiss} entries missing)')
                    fig, ax_hist = plt.subplots(
                        figsize = (6,4)
                    )
                    ax_hist.grid(True)
                    sns.histplot(data = df, x = v,
                                    kde=True, color = "green",
                                    ax = ax_hist)
                    fig.set_tight_layout(True)
                    st.write(fig)
                jumptwice()


        # Once all QC is done, we are are ready to finalise the process
        if st.button("Finished?"):
            #if not (clinical_conf & ph_conf & sex_conf & race_conf & fh_conf & rg_conf):
            if not (ph_conf & sex_conf & race_conf & fh_conf & rg_conf):
                st.error('Did you forget to confirm any of the steps above?')
                st.error("Please, tick all the boxes on the previous steps if the QC to meet GP2 standard format was successful")
            elif studycode == None:#st.session_state['keepcode'] == None:
                st.error('Please make sure you have selected the study code on the side bar')
                st.error('If you are unsure which study code is yours, please email us at cohort@gp2.org')
            else:
                # Update json file with master IDs
                #generategp2ids.update_masterids(ids_log, study_tracker)
                for idslog_tracker in st.session_state['all_ids']:
                   generategp2ids.update_masterids(idslog_tracker[0], idslog_tracker[1], studycode)

                email_ellie(studycode = studycode, activity = 'qc')
                    
                df = df.reset_index(drop=True)
                df['CustomOrder'] = df['sample_id'].apply(lambda x: original_order.index(x))
                df = df.sort_values(by='CustomOrder')
                df = df.drop(columns=['CustomOrder'])

                st.session_state['smqc'] = df
                
                st.markdown('<p class="medium-font"> CONGRATS, your sample manifest meets all the GP2 requirements. </p>', unsafe_allow_html=True )
                writeexcel = to_excel(df = df, 
                                      studycode = st.session_state['keepcode'], 
                                      mv = mani_vers,
                                      datatype = 'sm')
                st.download_button(label='📥 Download your QC sample manifest',
                                   data = writeexcel[0],
                                   file_name = writeexcel[1])

        jumptwice()
        st.markdown("<a href='#link_to_top'>Link to top</a>", unsafe_allow_html=True)