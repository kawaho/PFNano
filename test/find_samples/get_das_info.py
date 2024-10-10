import subprocess
import json

def get_das_info(query):
    '''Interface with das.py to get the query output.
    Could be done better, but this is time effective.
    Unfortunately the QL is more complicated than the 
    DBS one. '''
    
    das_command = [ 
        'dasgoclient',
        '--query=%s' % query,
        '--limit=0' 
        ]   
    p = subprocess.Popen(
        das_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
        )   
    out, err = p.communicate()
    das_exitcode = p.wait()
    
    if das_exitcode != 0:
        #sometimes das sends the crash message to stdout
        raise RuntimeError(f'crashed with error:\n{err}')
    return [str(i.strip(), 'utf-8') for i in out.split(b'\n') if i.strip()]

campaigns = {'2016preVFP':'RunIISummer20UL16MiniAOD*v2-106X*_preVFP*', '2016postVFP':'RunIISummer20UL16MiniAOD*v2-106X*v17*', '2017':'RunIISummer20UL17MiniAOD*v2-106X*', '2018':'RunIISummer20UL18MiniAOD*v2-106X*'}
missing_samples = {'2016preVFP':[], '2016postVFP':[], '2017':[], '2018':[]}

with open("MiniAOD_MC.json", 'r') as f:
  MC_names = json.load(f)
  print(MC_names)

for year, campaign in campaigns.items():
  print(f"----------------------Checking MC samples DAG for year {year}----------------------")
  allsamples = {}
  for MC_shorthand, MC_name in MC_names.items():
    sample = get_das_info(f"/*{MC_name}*/{campaign}/*")
    if sample:
      allsamples[MC_shorthand]=list(sample)
      print(f"{MC_name} Found!!")
      for subsample in sample:
          if 'PUForMUOVal' in subsample or 'Pilot' in subsample:
            allsamples[MC_shorthand].remove(subsample) 
          elif 'ext' in subsample:
            allsamples.setdefault(MC_shorthand+'_ext',[]).append(subsample)
            allsamples[MC_shorthand].remove(subsample) 
      if len(allsamples[MC_shorthand]) > 1: 
        print(f"!!{MC_name} has duplicates!!")
    else:
      print(f"{MC_name} is Missing!!")
      missing_samples[year].append(MC_name)
  with open("MiniAODUL_%s_MC.json"%year, 'w') as f:
    json.dump(allsamples, f, indent=4, sort_keys=True)
with open("MissingSamples.json", 'w') as f:
  json.dump(missing_samples, f, indent=4, sort_keys=True)

campaigns = {'2016':'Run2016*UL2016_MiniAODv2-v*', '2017':'Run2017*-UL2017_MiniAODv2-v*', '2018':'Run2018*-UL2018_MiniAODv2-v*'}
dataNames = ['SingleMuon']
for year, campaign in campaigns.items():
  print(f"----------------------Checking DATA samples DAG for year {year}----------------------")
  allsamples = {}
  for dataname in dataNames:
    sample = get_das_info("/%s/%s/*"%(dataname,campaign))
    if sample:
      for run in sample:
        runName = run.split("/")[2].split("_")[0]
        subsample = dataname+'_'+runName
        if subsample in allsamples: 
          allsamples[subsample].append(run)
        else:
          allsamples[subsample]=[run]
        print(f"{subsample} Found!!")
  if year=='2016':
    allsamplespostVFP,  allsamplespreVFP = {}, {}
    for shorthand2016, names2016 in allsamples.items():     
      if '2016F-UL' in shorthand2016 or '2016G-' in shorthand2016 or '2016H-' in shorthand2016:
        allsamplespostVFP[shorthand2016]= names2016
      else: 
        allsamplespreVFP[shorthand2016] = names2016
 
    with open("MiniAODUL_2016preVFP_data.json", 'w') as f:
      json.dump(allsamplespreVFP, f, indent=4, sort_keys=True)
    with open("MiniAODUL_2016postVFP_data.json", 'w') as f:
      json.dump(allsamplespostVFP, f, indent=4, sort_keys=True)
  else:
    with open("MiniAODUL_%s_data.json"%year, 'w') as f:
      json.dump(allsamples, f, indent=4, sort_keys=True)


