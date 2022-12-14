try:
    import streamlit as st
    import streamlit.components.v1 as stc

    # File Processing Pkgs
    import pandas as pd
    import numpy as np
    import base64
    import datetime as dt
    import matplotlib.pyplot as plt
    import seaborn as sns

except Exception as e:
    print("Some modules are not installed {}".format(e))


def jumptwice():
    st.write("##")
    st.write("##")

@st.cache
def read_file(data_file):
    if data_file.type == "text/csv":
        df = pd.read_csv(data_file)
    elif data_file.type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet":
        df = pd.read_excel(data_file, sheet_name=0)
    return (df)

def get_table_download_link(df, filetype = "CSV"):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in tsv or csv format
    out: href string
    """
    today = dt.datetime.today()
    version = f'{today.year}{today.month}{today.day}'

    study_code = df.study.unique()[0]
    if filetype == "CSV":
        csv = df.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
        href = f'<a href="data:file/csv;base64,{b64}"  download="{study_code}_sample_manifest_selfQC_{version}.csv">Download csv file</a>'
    else:
        tsv = df.to_csv(index=False, sep="\t")
        b64 = base64.b64encode(tsv.encode()).decode()  # some strings <-> bytes conversions necessary here
        href = f'<a href="data:file/tsv;base64,{b64}"  download="{study_code}_sample_manifest_selfQC_{version}.tsv">Download tsv file</a>'
    return href

def app():
    
    st.markdown("""
    <div id='linkto_top'></div>
    <style>
    .big-font {
        font-family:Helvetica; color:#0f557a; font-size:48px !important;
    }
    .medium-font {
        font-family:Arial; color:000000; font-size:18px;
    }
    </style>
    """, unsafe_allow_html=True)

    menu = ["For Fulgent", "For NIH"]
    choice = st.sidebar.selectbox("Menu",menu)
    
    ph_conf=''
    sex_conf=''
    race_conf = ''
    fh_conf=''
    
    cols = ['study', 'sample_id', 'sample_type',
            'DNA_volume', 'DNA_conc', 'r260_280',
            'Plate_name', 'Plate_position', 'clinical_id', 
            'study_arm', 'sex', 'race', 
            'age', 'age_of_onset', 'age_at_diagnosis', 'age_at_death', 'family_history',
            'region', 'comment', 'alternative_id1', 'alternative_id2']
    required_cols = ['study', 'sample_id', 'sample_type', 'clinical_id','study_arm', 'sex']
    allowed_samples = ['Blood (EDTA)', 'Blood (ACD)', 'Blood', 'DNA',
                        'DNA from blood', 'DNA from FFPE', 'RNA', 'Saliva',
                        'Buccal Swab', 'T-25 Flasks (Amniotic)', 'FFPE Slide',
                        'FFPE Block', 'Fresh tissue', 'Frozen tissue',
                        'Bone Marrow Aspirate', 'Whole BMA', 'CD3+ BMA', 'Other']
    fulgent_cols = ['DNA_volume', 'DNA_conc', 'Plate_name', 'Plate_position']
    
    data_file = st.sidebar.file_uploader("Upload Sample Manifest (CSV/XLSX)", type=['csv', 'xlsx'])

    file_type = ["CSV", "TSV"]
    output_choice = st.sidebar.selectbox("Select the output format of your sample manifest", file_type)
    
    # Process allowed_samples
    allowed_samples_strp = [samptype.strip().replace(" ", "") for samptype in allowed_samples]

    st.markdown('<p class="big-font">GP2 sample manifest self-QC</p>', unsafe_allow_html=True)
    #st.title("GP2 sample manifest self-QC web app")
    st.markdown('<p class="medium-font"> This is a web app to self-check the sample manifest </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Download the template from the link below. Once you open the link, go to "File"> "Download" as xlsx/csv). </p>', unsafe_allow_html=True )
    st.write('https://docs.google.com/spreadsheets/d/1ThpUVBCRaPdDSiQiKZVwFpwWbW8mZxhqurv7-jBPrM4')
    st.markdown('<p class="medium-font"> Please refer to the second tab (Dictionary) for instructions </p>', unsafe_allow_html=True)
    st.markdown('<p class="medium-font"> Once all the GP2 required columns are present in your manifest, please upload the sample manifest on the side bar for self QC.  </p>', unsafe_allow_html=True)
    
    #st.text('This is a web app to self-check the sample manifest')
    #st.text('The template download from the below link (go to "File"> "Download" as xlsx/csv)')
    
    #st.text('Please refer to the second tab (Dictionary) for instruction')
    st.text('')
    if data_file is not None:
        st.header("Data Check and self-QC")
        
        # read a file
        df = read_file(data_file)
        df['Genotyping_site'] = choice.replace('For ', '')
        if choice=='For Fulgent':
            required_cols = required_cols + fulgent_cols
        
        # missing col check
        missing_cols = np.setdiff1d(cols, df.columns)
        if len(missing_cols)>0:
            st.error(f'{missing_cols} are missing. Please use the template sheet')
            st.stop()
        else:
            st.text('Check column names--> OK')
            df_non_miss_check = df[required_cols].copy()

        # sample type check
        st.text('sample_type check')
        st.write(df.sample_type.astype('str').value_counts())
        not_allowed = np.setdiff1d(df.sample_type.unique(), allowed_samples)
        if len(not_allowed)>0:
            sampletype = df.sample_type
            sampletype = sampletype.str.strip().replace(" ", "")
            not_allowed_v2 = np.setdiff1d(sampletype.unique(), allowed_samples_strp)
            
            if len(not_allowed_v2) < len(not_allowed):
                st.text('We have found some undesired whitespaces in some sample type values.')
                    
                if len(not_allowed_v2)>0:
                    st.text('In addition, we have found unknown sample types')
                    st.error(f'sample_type: {not_allowed} not allowed.')
                    sample_list = '\n * '.join(allowed_samples)
                    st.text(f'Allowed sample list - \n * {sample_list}')
                    st.stop()

                else: # Map the stripped sample types back to the harmonised names
                    st.text('Processing whitespaces found in certain sample_type entries')
                    stype_map = dict(zip(allowed_samples_strp, allowed_samples))
                    newsampletype = sampletype.replace(stype_map)
                    df['sample_type'] = newsampletype
                    st.text('sample type after removing undesired whitespaces')
                    st.write(df.sample_type.astype('str').value_counts())

        
        # required columns checks
        if df_non_miss_check.isna().sum().sum()>0:
            st.error('There are some missing entries in the required columns. Please fill the missing cells ')
            st.text('First ~30 columns with missing data in any required fields')
            st.write(df_non_miss_check[df_non_miss_check.isna().sum(1)>0].head(20))
            st.stop()
        else:
            st.text('Check missing data in the required fields --> OK')
            sample_id_dup = df.sample_id[df.sample_id.duplicated()].unique()    
        
        # sample dup check
        if len(sample_id_dup)>0:
            st.text(f'Duplicated sample_id:{sample_id_dup}')
            st.error(f'Unique sample IDs are required (clinical IDs can be duplicated if replicated)')
            st.stop()
        else:
            st.text(f'Check sample_id duplicaiton --> OK')
            st.text(f'N of sample_id (entries):{df.shape[0]}')
            st.text(f'N of unique clinical_id : {len(df.clinical_id.unique())}')

        # # sample type check
        # st.text('sample_type check')
        # st.write(df.sample_type.astype('str').value_counts())
        # not_allowed = np.setdiff1d(df.sample_type.unique(), allowed_samples)
        # if len(not_allowed)>0:
        #     sampletype = df.sample_type
        #     sampletype = sampletype.str.strip().replace(" ", "")
        #     not_allowed_v2 = np.setdiff1d(sampletype.unique(), allowed_samples_strp)
            
        #     if len(not_allowed_v2) < len(not_allowed):
        #         st.text('We have found some undesired whitespaces in some sample type values.')
                    
        #         if len(not_allowed_v2)>0:
        #             st.text('In addition, we have found unknown sample types')
        #             st.error(f'sample_type: {not_allowed} not allowed.')
        #             sample_list = '\n * '.join(allowed_samples)
        #             st.text(f'Allowed sample list - \n * {sample_list}')
        #             st.stop()

        #         else: # Map the stripped sample types back to the harmonised names
        #             st.text('Processing whitespaces found in certain sample_type entries')
        #             stype_map = dict(zip(allowed_samples_strp, allowed_samples))
        #             newsampletype = sampletype.replace(stype_map)
        #             df['sample_type'] = newsampletype
        #             st.text('sample type after removing undesired whitespaces')
        #             st.write(df.sample_type.astype('str').value_counts())

        # study_arm --> Phenotype
        jumptwice()
        st.subheader('Create "Phenotype"')
        st.text('Count per study_arm')
        st.write(df.study_arm.astype('str').value_counts())
        arms=df.study_arm.dropna().unique()
        n_arms = st.columns(len(arms))
        phenotypes={}
        for i, x in enumerate(n_arms):
            with x:
                arm = arms[i]
                phenotypes[arm]=x.selectbox(f"[{arm}]: For QC, please pick the closest Phenotype",["PD", "Control", "Prodromal", "Other", "Not Reported"], key=i)
        df['Phenotype'] = df.study_arm.map(phenotypes)

        # cross-tabulation of study_arm and Phenotype
        st.text('=== Phenotype x study_arm===')
        xtab = df.pivot_table(index='Phenotype', columns='study_arm', margins=True,
                                values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)
        
        ph_conf = st.checkbox('Confirm Phenotype?')
        if ph_conf:
            st.info('Thank you')
        

        # sex for qc
        jumptwice()
        st.subheader('Create "sex_for_qc"')
        st.text('Count per sex group')
        st.write(df.sex.astype('str').value_counts())
        sexes=df.sex.dropna().unique()
        n_sexes = st.columns(len(sexes))
        mapdic={}
        for i, x in enumerate(n_sexes):
            with x:
                sex = sexes[i]
                mapdic[sex]=x.selectbox(f"[{sex}]: For QC, please pick a word below", 
                                    ["Male", "Female", "Intersex", "Unknown", "Other", "Not Reported"], key=i)
        df['sex_for_qc'] = df.sex.replace(mapdic)

        # cross-tabulation of study_arm and Phenotype
        st.text('=== sex_for_qc x sex ===')
        xtab = df.pivot_table(index='sex_for_qc', columns='sex', margins=True,
                                values='sample_id', aggfunc='count', fill_value=0)            
        st.write(xtab)
        
        sex_conf = st.checkbox('Confirm sex_for_qc?')
        if sex_conf:
            st.info('Thank you')

        # race for qc
        jumptwice()
        st.subheader('Create "race_for_qc"')
        st.text('Count per race (Not Reported = missing)')
        df['race_for_qc'] = df.race.fillna('Not Reported')
        st.write(df.race_for_qc.astype('str').value_counts())
        races = df.race.dropna().unique()
        nmiss = sum(pd.isna(df.race))

        if nmiss>0:
            st.text(f'{nmiss} entries missing race...')
            
        
        mapdic = {'Not Reported':'Not Reported'}
        for race in races:
            mapdic[race]=st.selectbox(f"[{race}]: For QC purppose, select the best match from the followings",
            ["American Indian or Alaska Native", "Asian", "White", "Black or African American", 
            "Multi-racial", "Native Hawaiian or Other Pacific Islander", "Other", "Unknown", "Not Reported"])
        df['race_for_qc'] = df.race_for_qc.map(mapdic)
        
        # cross-tabulation
        st.text('=== race_for_qc X race ===')
        dft = df.copy()
        dft['race'] = dft.race.fillna('_Missing')
        xtab = dft.pivot_table(index='race_for_qc', columns='race', margins=True,
                                values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)
        
        race_conf = st.checkbox('Confirm race_for_qc?')
        if race_conf:
            st.info('Thank you')
        

        # family history for qc
        jumptwice()
        st.subheader('Create "family_history_for_qc"')
        st.text('Count per family_history category (Not Reported = missing)')
        df['family_history_for_qc'] = df.family_history.fillna('Not Reported')
        st.write(df.family_history_for_qc.astype('str').value_counts())
        family_historys = df.family_history.dropna().unique()
        nmiss = sum(pd.isna(df.family_history))

        if nmiss>0:
            st.text(f'{nmiss} entries missing family_history')
        mapdic = {'Not Reported':'Not Reported'}

        if len(family_historys)>0:
            n_fhs = st.columns(len(family_historys))
            for i, x in enumerate(n_fhs):
                with x:
                    fh = family_historys[i]
                    mapdic[fh]=x.selectbox(f'[{fh}]: For QC, any family history?',['Yes', 'No', 'Not Reported'], key=i)
        df['family_history_for_qc'] = df.family_history_for_qc.map(mapdic).fillna('Not Assigned')

        # cross-tabulation 
        st.text('=== family_history_for_qc X family_history ===')
        dft = df.copy()
        dft['family_history'] = dft.family_history.fillna('_Missing')
        xtab = dft.pivot_table(index='family_history_for_qc', columns='family_history', margins=True,
                                values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)

        fh_conf = st.checkbox('Confirm family_history_for_qc?')
        if fh_conf:
            st.info('Thank you')

        # Region for qc
        jumptwice()
        st.subheader('Create "region_for_qc"')
        st.text('Count per region (Not Reported = missing)')
        df['region_for_qc'] = df.region.fillna('Not Reported')
        st.write(df.region_for_qc.astype('str').value_counts())
        regions = df.region.dropna().unique()
        nmiss = sum(pd.isna(df.region))
        if nmiss>0:
            st.text(f'{nmiss} entries missing for region')
        
        mapdic = {'Not Reported':'Not Reported'}

        if len(regions)>0:
            st.text('if ISO 3166-3 is available for the region, please provide')
            st.write('https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3')
            n_rgs = st.columns(len(regions))
            for i, x in enumerate(n_rgs):
                with x:
                    region = regions[i]
                    region_to_map = x.text_input(f'[{region}] in 3 LETTER (or NA)')
                    if len(region_to_map)>1:
                        mapdic[region]=region_to_map
        df['region_for_qc'] = df.region_for_qc.map(mapdic).fillna('Not Assigned')

        # cross-tabulation 
        st.text('=== region_for_qc X region ===')
        dft = df.copy()
        dft['region'] = dft.region.fillna('_Missing')
        xtab = dft.pivot_table(index='region_for_qc', columns='region', margins=True,
                                values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)

        rg_conf = st.checkbox('Confirm region_for_qc?')
        if rg_conf:
            if len(regions)!=(len(mapdic)-1):
                st.error('region not assigned.')
                st.text("Please, assing the region, and re upload the sample manifest")
                st.stop()
            else:
                st.info('Thank you')


        # Plate Info
        jumptwice()
        st.subheader('Plate Info')
        dft = df.copy()
        dft['Plate_name'] = dft.Plate_name.fillna('_Missing')
        xtab = dft.pivot_table(index='Plate_name', 
                            columns='study_arm', margins=True,
                            values='sample_id', aggfunc='count', fill_value=0)
        st.write(xtab)

        for plate in dft.Plate_name.unique():
            df_plate = dft[dft.Plate_name==plate].copy()
            df_plate_pos = df_plate.Plate_position
            # duplicated position check
            if plate!='_Missing':
                if len(df_plate_pos)>96:
                    st.error('Please make sure, N of samples on plate [{plate}] is =<96')
                    st.stop()
                dup_pos = df_plate_pos[df_plate_pos.duplicated()].unique()
                if len(dup_pos)>0:
                    st.error(f' !!!SERIOUS ERROR!!!  Plate position duplicated position {dup_pos} on plate [{plate}]')
                    st.stop()

        # Numeric values
        jumptwice()
        st.subheader('Numeric Values')
        numerics_cols = ['DNA_volume', 'DNA_conc', 'r260_280','age', 'age_of_onset', 'age_at_diagnosis', 'age_at_death']
        for v in numerics_cols:
            if df.dtypes[v] not in ['float64', 'int64']:
                st.error(f'{v} is not numeric')
                st.text("Please, make sure expected numeric columns are stored on a numeric format")
                st.text(f'Expected columns in numeric format - \n * {numerics_cols}')
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
                    #fig, ax = plt.subplots(figsize=(11,6))
                    #ax.hist(df[v].values, bins=20)
                    #st.pyplot(fig)
                    fig, ax_hist = plt.subplots(
                        figsize = (11,6)
                    )
                    ax_hist.grid(True)
                    sns.histplot(data = df, x = v,
                                    kde=True, color = "green",
                                    ax = ax_hist)
                    fig.set_tight_layout(True)
                    st.pyplot(fig)
                    #import plotly_express as px
                    #plot = px.histogram(x=v, data_frame=df)
                    #st.plotly_chart(plot)
                jumptwice()

        if st.button("Finished?"):
            if not (ph_conf & sex_conf & race_conf & fh_conf & rg_conf):
                st.error('Did you forget to confirm any of the steps above?')
                st.text("Please, tick all the boxes on the previous steps if the QC to meet GP2 standard format was successful")
            else:
                st.markdown('<p class="medium-font"> CONGRATS, your sample manifest meets all the GP2 requirements. </p>', unsafe_allow_html=True )
                st.markdown('<p class="medium-font"> Please, download it from the link below </p>', unsafe_allow_html=True )
                if output_choice == "CSV":
                    st.markdown(get_table_download_link(df), unsafe_allow_html=True)
                else:
                    st.markdown(get_table_download_link(df, filetype=output_choice), unsafe_allow_html=True)
                jumptwice()
                st.markdown('<p class="medium-font"> Now, click the link below to go back to the top of the page and change to the Upload tab to deposit the QC sample manifest on the GP2 storage system</p>', unsafe_allow_html=True )
                st.markdown("<a href='#linkto_top'>Link to top</a>", unsafe_allow_html=True)