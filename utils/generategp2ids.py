import json
import ijson
import numpy as np
import pandas as pd
from google.cloud import storage
import streamlit as st
import datetime as dt

#@st.cache(hash_funcs={'_json.Scanner': hash})
#@st.experimental_memo()




def update_masterids(ids_log, ids_dict, scode):
    client = storage.Client()
    bucket = client.get_bucket('eu-samplemanifest')
    blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER.json')
    #blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER_for_monogenic.json')
    #blob = bucket.blob('IDSTRACKER/TESTS_GP2IDSMAPPER.json')
    
    today = dt.datetime.today()
    version = f'{today.year}{today.month}{today.day}'
    archive_blob = bucket.blob("IDSTRACKER/ARCHIVE/{v}_{s}_GP2IDSMAPPER.json".format(v = version, s = scode))

    # Create security copy
    with blob.open("r") as fp:
        masterids = json.load(fp)
    with archive_blob.open("w") as fp:
        json.dump(masterids, fp)
    if bool(ids_dict):
        for study, newdata in ids_log.items():
            masterids[study].update(newdata)
    else:
        for study, newdata in ids_log.items():
            masterids[study] = newdata
    
    # Update MASTERJSON with new IDs
    with blob.open("w") as fp:
        json.dump(masterids, fp)
    
    return(masterids)

#@st.cache
def master_key(studies):
    # ACCESS MASTERGP2IDS_JSON IN GP2 BUCKET
    client = storage.Client()
    bucket = client.get_bucket('eu-samplemanifest')
    blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER.json')
    #blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER_for_monogenic.json')
    #blob = bucket.blob('IDSTRACKER/TESTS_GP2IDSMAPPER.json')

    ids_tracker = {}
    for cohort in studies:
        with blob.open("r") as f:
        #for cohort in studies:
            for record in ijson.items(f, cohort):
                ids_tracker.update({cohort : record})
    
    return(ids_tracker)


def master_keyv2(studies):
    # ACCESS MASTERGP2IDS_JSON IN GP2 BUCKET
    client = storage.Client()
    bucket = client.get_bucket('eu-samplemanifest')
    blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER.json')
    #blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER_for_monogenic.json')
    #blob = bucket.blob('IDSTRACKER/TESTS_GP2IDSMAPPER.json')
    ids_tracker = {}
    with blob.open("r") as f:
        for k, v in ijson.kvitems(f, ''):
            if k in studies:
                ids_tracker.update({k:v})
    return(ids_tracker)


def master_remove(studies, data):
    client = storage.Client()
    bucket = client.get_bucket('eu-samplemanifest')
    #blob = bucket.blob('IDSTRACKER/TESTS_GP2IDSMAPPER.json')
    #blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER_for_monogenic.json')
    blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER.json')
    
    with blob.open("r") as fp:
        masterids = json.load(fp)
    
    masterids_copy = masterids.copy()
    for study in studies:
        study_master = masterids_copy.pop(study)
        if len(data['sample_id'].to_list()) == len(list(study_master.keys())):
            removestudy = masterids.pop(study)
        else:
            ids_remove = list(set(data['sample_id'].to_list()) & set(list(study_master.keys())))
            for k in ids_remove:
                removeid = study_master.pop(k, None)
    
    with blob.open("w") as fp:
        json.dump(masterids, fp) # Update the json file
    
    return(masterids)


def getgp2idsv2(dfproc, n, study_code):
    dfproc = dfproc.sort_values('sample_id')
    df_dups = dfproc[dfproc.duplicated(keep=False, subset=['clinical_id'])].sort_values('clinical_id').reset_index(drop = True).copy()
    
    if df_dups.shape[0]>0:
        dupids_mapper = dict(zip(df_dups.clinical_id.unique(),
                            [num+n for num in range(len(df_dups.clinical_id.unique()))]))
        
        df_dup_chunks = []
        for clin_id, gp2id in dupids_mapper.items():
            df_dups_subset = df_dups[df_dups.clinical_id==clin_id].copy()
            df_dups_subset['GP2ID'] = [f'{study_code}_{gp2id:06}' for i in range(df_dups_subset.shape[0])]
            df_dups_subset['SampleRepNo'] = ['s'+str(i+1) for i in range(df_dups_subset.shape[0])]
            df_dups_subset['GP2sampleID'] = df_dups_subset['GP2ID'] + '_' + df_dups_subset['SampleRepNo']
            df_dup_chunks.append(df_dups_subset)
        df_dups_wids = pd.concat(df_dup_chunks)

    df_nodups = dfproc[~dfproc.duplicated(keep=False, subset=['clinical_id'])].sort_values('clinical_id').reset_index(drop = True).copy()

    if not df_nodups.empty():
        if df_dups.shape[0]>0:
            n =  len(list(dupids_mapper.values())) + n
        else:
            n = n

        uids = [str(id) for id in df_nodups['sample_id'].unique()]
        mapid = {}
        for uid in uids:
            mapid[uid]= n
            n += 1
        df_nodups_wids = df_nodups.copy()
        df_nodups_wids['uid_idx'] = df_nodups_wids['sample_id'].map(mapid)
        df_nodups_wids['GP2ID'] = [f'{study_code}_{i:06}' for i in df_nodups_wids.uid_idx]
        df_nodups_wids['uid_idx_cumcount'] = df_nodups_wids.groupby('GP2ID').cumcount() + 1
        df_nodups_wids['GP2sampleID'] = df_nodups_wids.GP2ID + '_s' + df_nodups_wids.uid_idx_cumcount.astype('str')
        df_nodups_wids['SampleRepNo'] = 's' + df_nodups_wids.uid_idx_cumcount.astype('str')
        df_nodups_wids.drop(['uid_idx','uid_idx_cumcount'], axis = 1, inplace = True)

        if df_dups.shape[0]>0:
            df_newids = pd.concat([df_dups_wids, df_nodups_wids])
        else:
            df_newids = df_nodups_wids
    else:
        df_newids = df_dups_wids
    
    return(df_newids)

def assign_unique_gp2clinicalids(df, clinicalid_subset):

    if isinstance(clinicalid_subset, pd.Series):
        clinicalid_subset = clinicalid_subset.to_frame().T

    sampleid = clinicalid_subset.sort_values(by=['master_GP2sampleID'])\
                                .reset_index(drop = True)\
                                .dropna(subset=['master_GP2sampleID'], axis = 0)
    sampleid = sampleid.loc[sampleid.index[-1], 'master_GP2sampleID'].split("_")
    getuniqueid = sampleid[0] + "_" + sampleid[1]
    get_sidrepno = int(sampleid[2].replace("s","")) + 1

    index_modify = clinicalid_subset['index'].unique() #clinicalid_subset[clinicalid_subset['GP2sampleID'].isnull()] #.index
    assign_gp2sampleid = [getuniqueid + "_s" + str(get_sidrepno + i) for i in range(len(index_modify))]
    df.loc[index_modify, 'GP2sampleID'] = assign_gp2sampleid
    getnewidrows = df.loc[index_modify].copy()
    return (getnewidrows)



# def update_masterids(ids_log, ids_dict):
#     client = storage.Client()
#     bucket = client.get_bucket('eu-samplemanifest')
#     #blob = bucket.blob('IDSTRACKER/GP2IDSMAPPER.json')
#     blob = bucket.blob('IDSTRACKER/TESTS_GP2IDSMAPPER.json')
#     with blob.open("r") as fp:
#         masterids = json.load(fp)
#     if bool(ids_dict):
#         for study, newdata in ids_log.items():
#             masterids[study].update(newdata)
#     else:
#         for study, newdata in ids_log.items():
#             masterids[study] = newdata
#     with blob.open("w") as fp:
#         json.dump(masterids, fp)
#     return(masterids)


# def getgp2ids(df_updateids, n, uids, study_code):
#     mapid = {}
#     for uid in uids:
#         mapid[uid]= n
#         n += 1
#     mydf = df_updateids.copy()
#     mydf['uid_idx'] = mydf['sample_id'].map(mapid)
#     mydf['GP2ID'] = [f'{study_code}_{i:06}' for i in mydf.uid_idx]
#     mydf['uid_idx_cumcount'] = mydf.groupby('GP2ID').cumcount() + 1
#     mydf['GP2sampleID'] = mydf.GP2ID + '_s' + mydf.uid_idx_cumcount.astype('str')
#     mydf['SampleRepNo'] = 's' + mydf.uid_idx_cumcount.astype('str')
#     mydf.drop(['uid_idx','uid_idx_cumcount'], axis = 1, inplace = True)
#     return(mydf)


# @st.cache
# def blobbucket(bname = 'eu-samplemanifest'):
#     client = storage.Client()
#     bucket = client.get_bucket(bname)
#     blob = bucket.blob('IDSTRACKER/CLINICALGP2IDS_MAPPER_20230718.json')
#     return blob


# @st.cache
# def update_masterids(blob, ids_log, ids_dict):
#     #client = storage.Client()
#     #bucket = client.get_bucket('eu-samplemanifest')
#     #blob = bucket.blob('IDSTRACKER/CLINICALGP2IDS_MAPPER_20230718.json')
#     with blob.open("r") as fp:
#         masterids = json.load(fp)
    
#     #if 'ids_tracker' in globals():
#     if bool(ids_dict):
#         for study, newdata in ids_log.items():
#             masterids[study].update(newdata)
#     else:
#         for study, newdata in ids_log.items():
#             masterids[study] = newdata
    
#     with blob.open("w") as fp:
#         json.dump(masterids, fp)
    
#     print("MASTER IDS UPDATED")


# @st.cache
# def master_key(blob, studies):
#     # ACCESS MASTERGP2IDS_JSON IN GP2 BUCKET
#     #client = storage.Client()
#     #bucket = client.get_bucket('eu-samplemanifest')
#     #blob = bucket.blob('IDSTRACKER/CLINICALGP2IDS_MAPPER_20230718.json')
#     ids_tracker = {}
#     with blob.open("r") as f:
#         for cohort in studies:
#             st.write("IN")
#             for record in ijson.items(f, cohort):
#                 ids_tracker.update({cohort : record})
#     return(ids_tracker)