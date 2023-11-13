
geo_list2=tibble::tribble(
           ~GEO_id,                    ~Disease,    ~PathID, ~Normal, ~Case,   ~Pubmed,                              ~Tissue,
        "GSE14762",      "Renal cell carcinoma", "hsa05211",     12L,    9L, 19252501L,                             "Kidney",
         "GSE6344",      "Renal cell carcinoma", "hsa05211",     10L,   10L, 17699851L,                     "Clear cell RCC",
         "GSE1297",       "Alzheimer's disease", "hsa05010",      9L,    7L, 14769913L,                    "Hippocampal CA1",
       "GSE5281EC",       "Alzheimer's disease", "hsa05010",     13L,   10L, 17077275L,           "Brain, Entorhinal Cortex",
      "GSE5281HIP",       "Alzheimer's disease", "hsa05010",     13L,   10L, 17077275L,                 "Brain, hippocampus",
      "GSE5281VCX",       "Alzheimer's disease", "hsa05010",     12L,   19L, 17077275L,       "Brain, primary visual cortex",
        "GSE65144",            "Thyroid cancer", "hsa05216",     13L,   12L, 25675381L,                            "Thyroid",
        "GSE58545",            "Thyroid cancer", "hsa05216",     18L,   27L, 26625260L,                            "Thyroid",
        "GSE41011",         "Colorectal cancer", "hsa05210",     12L,   19L,        NA,                              "Colon",
        "GSE55945",           "Prostate cancer", "hsa05215",      7L,   12L, 19737960L,                           "Prostate",
        "GSE26910",           "Prostate cancer", "hsa05215",      6L,    6L, 21611158L,                           "Prostate",
         "GSE8762",      "Huntington's disease", "hsa05016",     10L,   12L, 17724341L,                         "Lymphocyte",
        "GSE24250",      "Huntington's disease", "hsa05016",      6L,    8L, 21969577L,        "Venous cellular whole blood",
        "GSE73655",      "Huntington's disease", "hsa05016",      7L,   13L, 26756592L,               "Subcutaneous adipose",
        "GSE37517",      "Huntington's disease", "hsa05016",      5L,    8L, 22748968L,                   "Neural stem cell",
    "GSE14924_CD4",    "Acute Myeloid Leukemia", "hsa05221",     10L,   10L, 19710498L,                         "CD4 T Cell",
        "GSE15471",         "Pancreatic cancer", "hsa05212",     35L,   35L, 19260470L,                           "Pancreas",
        "GSE16515",         "Pancreatic cancer", "hsa05212",     15L,   15L, 19732725L,                           "Pancreas",
        "GSE28735",         "Pancreatic cancer", "hsa05212",     45L,   45L, 23918603L,                           "Pancreas",
        "GSE20153",       "Parkinson's disease", "hsa05012",      8L,    8L, 20926834L, "Blymphocytes from peripheral blood",
        "GSE19587",      "Parkinson's disease", "hsa05012",     10L,   12L, 20837543L,                              "Brain",
        "GSE26887", "Type II diabetes mellitus", "hsa04930",      5L,    7L, 22427379L,                     "Left ventricle",
         "GSE7305",        "Endometrial cancer", "hsa05213",     10L,   10L, 17640886L,         "Endometrium/Ovarian tissue",
        "GSE36389",        "Endometrial cancer", "hsa05213",      7L,   13L,        NA,                        "Endometrium"
    )%>%as.data.frame()


geo_list3=tibble::tribble(
           ~GEO, ~KO_gene, ~Impacted_path_number, ~Normal, ~Case,   ~Pubmed,                                           ~Tissue,
     "GSE22873",  "Myd88",                24L,     11L,    8L, 22075646L,                                              "Liver",
    "GSE70302a",   "Il1a",                21L,      4L,    4L, 26224856L,                                        "Spinal cord",
    "GSE70302b",   "Il1b",                41L,      4L,    4L, 26224856L,                                        "Spinal cord",
     "GSE58120",    "Il2",                20L,      6L,    6L, 25652593L,                            "Myeloid dendritic cells",
     "GSE46211", "Tgfbr2",                23L,     12L,    6L, 24496627L,                            "Anterior palatal tissue",
    "GSE138957",   "Akt1",                97L,      3L,    3L, 32937140L,                                         "Cell lines",
     "GSE85754",  "Fgfr1",                16L,      4L,    4L, 28433771L,                                             "Breast",
     "GSE88799", "Crebbp",                28L,      4L,    5L, 28069569L,                                             "B cell",
      "GSE4451",    "Met",                21L,      6L,    6L, 16710476L,                                              "Liver"
    )%>%as.data.frame()

# "Bhlhe40", "Id3", "Dusp5", "Onecut1", "Neurod1"
pos_res=lapply(c("Myd88",  "Il1a", "Il1b", "Il2", "Tgfbr2","Akt1","Fgfr1","Crebbp","Met"),
       \(i)filter(mmu_kegg_pathway$all_org_gene,gene_symbol==i)%>%select(pathway_id)%>%mutate(mutate=i))%>%
    do.call(rbind,.)

get_sig_from_all=function(comparison_res_m){
    enrich_sig_m=list(
        #RS_m=filter(comparison_res_m$RS_m,(ReporterScore)>1.64,p.adjust<0.01)%>%pull(ID),
        GRSA=filter(comparison_res_m$RS_d,abs(ReporterScore)>1.64,p.value<0.05)%>%pull(ID),
        Fisher=filter(comparison_res_m$fisher_res,p.value<0.05)%>%pull(ID),
        CP=filter(comparison_res_m$enrich_res,p.value<0.05)%>%pull(ID),
        GSEA=filter(comparison_res_m$gsea_res@result,pvalue<0.05)%>%pull(ID)
        #GSA=filter(comparison_res_m$gsa_res,p.value<0.05)%>%pull(ID)
    )
    enrich_sig_m
}

get_sig_from_all2=function(comparison_res_m){
    enrich_sig_m=list(
        #RS_m=filter(comparison_res_m$RS_m,(ReporterScore)>1.64,p.adjust<0.01)%>%pull(ID),
        GRSA=filter(comparison_res_m$RS_d,abs(ReporterScore)>1.64,p.adjust<0.05)%>%pull(ID),
        Fisher=filter(comparison_res_m$fisher_res,p.adjust<0.05)%>%pull(ID),
        CP=filter(comparison_res_m$enrich_res,p.adjust<0.05)%>%pull(ID),
        GSEA=filter(comparison_res_m$gsea_res@result,p.adjust<0.05)%>%pull(ID)
        #GSA=filter(comparison_res_m$gsa_res,p.adjust<0.05)%>%pull(ID)
    )
    enrich_sig_m
}

load_geo=function(id="GSE65144"){
    id2=gsub("(GSE[0-9]+).*", "\\1", id)
    geo_tmp=readRDS(paste0("~/Documents/R/ReporterScore/article/refer data/GEO_data/tabs/",id2,".RDS"))

    switch (id,
            "GSE65144" = {meta_tmp=geo_tmp$meta%>%select(`tissue type:ch1`)%>%rename("group"=1)},
            "GSE33126" = {meta_tmp=geo_tmp$meta%>%select(`t=tumor; n=normal:ch1`)%>%rename("group"=1)%>%
                transmute(group=ifelse(group=="T","Tumor","Normal"))},
            "GSE41011" = {meta_tmp=geo_tmp$meta%>%select(`tissue:ch1`)%>%rename("group"=1)},
            #新加的数据
            "GSE5281EC" = {meta_tmp=geo_tmp$meta%>%
                select(`characteristics_ch1.8`,`characteristics_ch1.1`)%>%rename("group"=1,"type"=2)%>%
                filter(grepl("EC",type))%>%mutate(group=gsub("disease state","Disease State",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE5281HIP" = {meta_tmp=geo_tmp$meta%>%
                select(`characteristics_ch1.8`,`characteristics_ch1.1`)%>%rename("group"=1,"type"=2)%>%
                filter(grepl("HIP",type))%>%mutate(group=gsub("disease state","Disease State",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE5281VCX" = {meta_tmp=geo_tmp$meta%>%
                select(`characteristics_ch1.8`,`characteristics_ch1.1`)%>%rename("group"=1,"type"=2)%>%
                filter(grepl("VCX",type))%>%mutate(group=gsub("disease state","Disease State",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE1297" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)%>%
                filter(group%in%c("group: Severe","group: Control"));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE14762" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)},
            "GSE6344" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1`)%>%rename("group"=1)%>%mutate(group=gsub("(\\w+), PT.*", "\\1",group))},
            "GSE58545" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1`)%>%rename("group"=1)%>%mutate(group=gsub("type of thyroid tissue: ", "",group))},
            "GSE55945" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%
                mutate(group=ifelse(grepl("Normal",group),"Normal","Prostate Cancer"))},
            "GSE26910" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%
                mutate(group=gsub(" \\d+","",group))%>%filter(group%in%c("prostate normal","prostate tumor"));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE8762" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1`)%>%rename("group"=1)%>%
                mutate(group=ifelse(grepl("control",group),"control","HD patient"))},
            "GSE24250" = {meta_tmp=geo_tmp$meta%>%select(`source_name_ch1`)%>%rename("group"=1)},
            "GSE73655" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%
                mutate(group=ifelse(grepl("control",group),"control","HD patient"))},
            "GSE37517" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)},
            "GSE14924_CD4" = {meta_tmp=geo_tmp$meta%>%
                select(`title`,`source_name_ch1`)%>%rename("group"=1,"type"=2)%>%
                filter(grepl("CD4",type))%>%mutate(group=gsub("_CD4_.*","",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE15471" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)},
            "GSE16515" = {meta_tmp=geo_tmp$meta%>%select(`source_name_ch1`)%>%rename("group"=1)%>%
                mutate(group=ifelse(grepl("Normal",group),"Pancreatic-Normal","Pancreatic-Tumor"))},
            "GSE28735" = {meta_tmp=geo_tmp$meta%>%select(`source_name_ch1`)%>%rename("group"=1)%>%mutate(group=gsub(", patient \\d+","",group))},
            "GSE20153" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)},
            "GSE19587" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)%>%mutate(group=gsub("subject: (\\w+) .*","\\1",group))},
            "GSE26887" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.3`)%>%rename("group"=1)%>%
                mutate(group=gsub("disease state: *","",group))%>%
                filter(group%in%c("CONTROL","DIABETIC, HEART FAILURE"));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE7305" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%
                mutate(group=gsub(" \\d+","",group))},
            "GSE36389" = {meta_tmp=geo_tmp$meta%>%select(`source_name_ch1`)%>%rename("group"=1)%>%mutate(group=gsub(" G\\d+","",group))},
            "GSE58120" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.2`)%>%rename("group"=1)%>%mutate(group=gsub("genotype: ","",group))},
            "GSE22873" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)%>%
                mutate(group=ifelse(grepl("Myd88 null",group),"Myd88 null","WT"))},
            "GSE70302a" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)%>%filter(!grepl("Il1b",group))%>%mutate(group=gsub("disease state: ","",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE70302b" = {meta_tmp=geo_tmp$meta%>%select(`characteristics_ch1.1`)%>%rename("group"=1)%>%filter(!grepl("Il1a",group))%>%mutate(group=gsub("disease state: ","",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE46211" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%mutate(group=ifelse(grepl("Control",group),"Control","Mutant"))},
            "GSE138957" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%mutate(group=ifelse(grepl("\\+",group),"Control","Mutant"))},
            "GSE85754" = {meta_tmp=geo_tmp$meta%>%select(`title`)%>%rename("group"=1)%>%mutate(group=ifelse(grepl("KO",group),"FGFR1 KO","WT"))},
            "GSE88799" = {meta_tmp=geo_tmp$meta%>%select(`source_name_ch1`)%>%rename("group"=1)%>%filter(!grepl("HET",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]},
            "GSE4451" = {meta_tmp=geo_tmp$meta%>%select(`source_name_ch2`)%>%rename("group"=1)%>%mutate(group=gsub(" - Cy[3,5]","",group))%>%filter(grepl("Untreated",group));
            geo_tmp$GSE_expr=geo_tmp$GSE_expr[,rownames(meta_tmp)]}

    )
    #View(geo_tmp$meta)
    geo_tmp$meta=meta_tmp
    geo_tmp
}

#10.29更新，并且KO_gsea用Z_score排序分析,比较富集方法，选取p-value最小，logFC最大的前10%的gene做富集
hsa_gene=unique(hsa_kegg_pathway$all_org_gene$gene_symbol)
mmu_gene=unique(mmu_kegg_pathway$all_org_gene$gene_symbol)

test_tests4=function(KO_abundance,group,metadata,name,method="wilcox.test",p.adjust=T){
    if(file.exists(paste0(name,"_methods_res_d.RDS"))) {
        methods_res_d=readRDS(paste0(name,"_methods_res_d.RDS"))
    }
    else {
        if(p.adjust){
            methods_res_d=reporter_score(KO_abundance,group = group,
                                         metadata = metadata,mode="directed",method = method,feature = "gene",type = "mmu",threads = 1)
        }
        else {
            methods_res_d=reporter_score(KO_abundance,group = group,p.adjust.method1 = "none",
                                         metadata = metadata,mode="directed",method = method,feature = "gene",type = "mmu",threads = 1)
        }
        saveRDS(methods_res_d,file = paste0(name,"_methods_res_d.RDS"))
    }

    #reporter_score
    #fisher
    res.dt=methods_res_d$ko_stat
    # if(attributes(methods_res_d$reporter_s)$type=="hsa")res.dt=filter(res.dt,KO_id%in%hsa_gene)
    # if(attributes(methods_res_d$reporter_s)$type=="mmu")res.dt=filter(res.dt,KO_id%in%mmu_gene)

    add_mini=NULL
    if(!"logFC"%in%colnames(res.dt)){
        message("No logFC in the data.frame, calculate.")
        vs_group=grep("average",colnames(res.dt),value = T)
        if(length(vs_group)!=2)stop("logFC only available for two groups")
        tmp=c(res.dt[,vs_group[1]],res.dt[,vs_group[2]])

        if (is.null(add_mini))
            add_mini = min(tmp[tmp > 0]) * 0.05
        res.dt$logFC=log2((res.dt[,vs_group[2]]+add_mini)/(res.dt[,vs_group[1]]+add_mini))
    }

    sig_gene=filter(res.dt,p.value<0.05)%>%top_n(400,abs(logFC))%>%pull(KO_id)
    res.dt=mutate(res.dt,origin_p.adjust=ifelse(KO_id%in%sig_gene,0.04,0.1))
    fisher_res=KO_fisher(res.dt,padj_threshold = 0.05,modulelist = methods_res_d$modulelist)
    #enricher
    enrich_res=KO_enrich(res.dt,padj_threshold = 0.05,modulelist = methods_res_d$modulelist)

    #GESA
    set.seed(1234)
    gsea_res=KO_gsea(methods_res_d,padj_threshold = 1.1,weight = "Z_score")
    #GSA
    gsa_res=KO_gsa(methods_res_d)

    comparison_res_m=c(list(RS_d=methods_res_d$reporter_s),list(fisher_res=fisher_res),
                       list(enrich_res=enrich_res),gsea_res=gsea_res,list(gsa_res=gsa_res))

    saveRDS(comparison_res_m,file = paste0(name,"_comparison5.RDS"))
}

similar=\(a,b){tmp=two_set(a,b);tmp[2]/sum(tmp)}
similar_each=function(ls){
    n=length(ls)
    sim_mat=matrix(0,nrow = n,ncol = n)
    rownames(sim_mat)=colnames(sim_mat)=names(ls)
    for (i in 1:n) {
        for (j in i:n) {
            sim_mat[i,j]= sim_mat[j,i]=similar(ls[[i]],ls[[j]])
        }
    }
    as.data.frame(sim_mat)
}

