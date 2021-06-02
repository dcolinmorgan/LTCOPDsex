j=-1
jeff=pd.DataFrame()
# titleArray=['Binary network 100 genes','Continous network 1 100 genes','Continous network 2 100 genes','Continous network 3 100 genes']
traces= glob.glob('analyses/LTRC/bench/*2021-06-01*.txt')
traces=np.sort(traces)
for j,file in enumerate(traces):


    # for net in #[KRCCNetBin1,KRCCNetCont1,KRCCNetCont2,KRCCNetCont3]:
    j=j+1
    if '.txt' in(file):
        net=pd.read_csv(file,sep='\t',usecols=[3,8,9])
    else:
        net=pd.read_csv(file,sep=',',usecols=[3,8,9])
    net=net.sort_values('adj.P.Val').dropna()[0:500]
    c=net['GencodeBasicV12_NAME'].dropna().astype(object)
    genes_str='\n'.join(c)
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    query_string = '?userListId=%s&backgroundType=%s'

    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)

    user_list_id = data['userListId']
    gene_set_library = 'GO_Biological_Process_2018'
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text);
    processBio=[]
    pvals=[]
    for i in range(1,10):
        processBio.append(data['GO_Biological_Process_2018'][i][1])
        pvals.append(data['GO_Biological_Process_2018'][i][2])
    plt.figure()
    y_pos = np.arange(len(processBio))
    plt.barh(y_pos, -np.log10(pvals))
    plt.yticks(ticks=y_pos,labels=processBio);
    plt.gca().invert_yaxis()  
    plt.xlabel('-log(p-value)')
    plt.title(file.split('/')[3])
    pd.DataFrame(processBio).to_csv('~/analyses/LTRC/bench/GEO/'+file.split('/')[3])
#     except ValueError:
#         pass

    jeff=jeff.append(['~/analyses/LTRC/bench/GEO/'+file.split('/')[3]])
    jeff=jeff.append(pd.DataFrame(processBio))
jeff.to_csv('~/analyses/LTRC/bench/overall_enrich.txt',sep='\t')
