import pandas as pd
from collections import Counter

def annotate_raw(annotation_df, loci, one_length_print=False):
    # if one_length_print is true, we don't care about situation when the annotation file has only one entries for a hit, if that hit is a gene, if we are doing the double checking step with the promoter extended annotation file

    new_df=pd.DataFrame() #output dataframe with corresponding region for each high Fst locus

    for index, row in loci.iterrows(): #go through all high Fst loci
        ch = row[0].split("_")[:2]
        ch = "_".join(ch) #get chromosome
        pos = int(row[0].split("_")[2]) #get position for loci
        temp=annotation_df[(annotation_df['chr']==ch) & (annotation_df['start']<=pos) & (annotation_df['stop']>=pos)] #get rows in annotation file where our loci is between start and stop
        #print(len(temp))
        if len(temp) == 0:
            #print("went into len == 0")
            #locus is not annotated
            new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "not_annot", "gene" : "NA"}, ignore_index = True)
        if len(temp) == 1:
            #print("went into len ==1")
            #the annotation file ususally has mulitple lines belonging to a certain annotated region
            #e.g. for a gene, exon information is in a separate line
            #this is rare so I print if this is the case and inspect individual cases
            if one_length_print:
                if list(temp['type'])[0] != "gene":
                    print("ERROR: Only one line in annotation file, type of annotation:")
                    print(list(temp['type'])[0])
            else:
                if list(temp['type'])[0] == "pseudogene":
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "pseudogene", "gene" : "NA"}, ignore_index = True)
                else:
                    print("ERROR: Only one line in annotation file, type of annotation:")
                    print(list(temp['type'])[0])

        if (len(temp) > 1 or one_length_print) and (len(temp) != 0):
            #print("went into main if, len > 1")
            #checking if it is a known region, and if yes, if there are overlapping regions
            ids=[]
            for i in range(len(temp)):
                t=temp["info"].iloc[i].split(";")
                res = list(filter(lambda x: "LOC" in x, t))
                if len(res) == 0:
                    res = list(filter(lambda x: "GeneID:" in x, t))
                    if len(res) != 0:
                        ID = res[0].split("GeneID:",1)[1]
                        ids.append(ID)
                else:
                    ID =res[0].split("LOC",1)[1]
                    ids.append(ID)
            #ids is a list of IDs we found in the rows where our loci is between start and stop
            uni=len(Counter(ids))
            if uni > 1:
                #if there is more than 1 kind of gene ID, there are overlapping genes, investigate cases separately
                new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "genes_overlap", "gene": "NA"}, ignore_index = True)

            gene_name = ids[0] #technically it might not be a gene

            if uni == 1:
                #if there is only one kind of region
                #checking what kind of region
                if temp['type'].str.contains('gene').any():
                    if temp['type'].str.contains('pseudogene').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "pseudogene", "gene":"NA" },ignore_index = True)
                    # snRNA, tRNA, lnc_RNA, snoRNA, rRNA, miRNA also have line with "gene"
                    elif temp['type'].str.contains('snRNA').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "snRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('tRNA').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "tRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('lnc_RNA').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "lnc_RNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('snoRNA').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "snoRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('rRNA').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "rRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('miRNA').any():
                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "miRNA", "gene":ids[0] },ignore_index = True)
                    else:
                        if temp['type'].str.contains('exon').any():
                            if temp['type'].str.contains('CDS').any():
                                #in exon
                                new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "exon", "gene":f"LOC{ids[0]}" },ignore_index = True)
                            else:
                                direction=temp['s2'].iloc[0]
                                e_start=temp[temp["type"]=="exon"].iloc[0].start
                                e_stop=temp[temp["type"]=="exon"].iloc[0].stop
                                g_start=temp[temp["type"]=="gene"].iloc[0].start
                                g_stop=temp[temp["type"]=="gene"].iloc[0].stop
                                if direction == "+":
                                    if e_start == g_start:
                                        #print('5prime1')
                                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "5'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                                    elif e_stop == g_stop:
                                        #print("3prime1")
                                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "3'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                                    else:
                                        #print("ERROR1") #alternative UTR, inspect separately
                                        #print(temp)
                                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "alternative UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                                elif direction == "-":
                                    if e_start == g_start:
                                        #print('3prime2')
                                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "3'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                                    elif e_stop == g_stop:
                                        #print("5prime2")
                                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "5'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                                    else:
                                        #print("ERROR2") #alternative UTR, inspect separately
                                        new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "alternative UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                                else:
                                    print("ERROR3") #this should never return
                        else:
                            #in intron
                            new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "intron", "gene":f"LOC{ids[0]}" },ignore_index = True)
                else:
                    print("not gene")
    return new_df

def annotate_raw_region(annotation_df, loci): # NOTE: didn't end up using it
#use this function instead if not using single SNPs, instead regions, and want to see if region is superset of gene in annotation file
  new_df=pd.DataFrame() #output dataframe with corresponding region for each high Fst locus

  for index, row in loci.iterrows(): #go through all high Fst loci
      ch = row[0]
      #print(ch)
      pos1 = int(row[3]) #get start of region
      #print(pos1)
      pos2 = int(row[4]) #get end of region
      pos=str(pos1)+"_"+str(pos2)
      temp=annotation_df[(annotation_df['chr']==ch) & (pos1<=annotation_df['start']) & (pos2>=annotation_df['stop'])] #get rows in annotation file where our loci is between start and stop
      if len(temp) == 0:
          #locus is not annotated
          new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "not_annot", "gene" : "NA"},
                          ignore_index = True)
      if len(temp) == 1:
          #the annotation file ususally has mulitple lines belonging to a certain annotated region
          #e.g. for a gene, exon information is in a separate line
          #this is rare so I print if this is the case and inspect individual cases
          print("ERROR: Only one line in annotation file, type of annotation:")
          print(temp['type'])

      if len(temp) > 1:
          #checking if it is a known region, and if yes, if there are overlapping regions
          ids=[]
          for i in range(len(temp)):
              t=temp["info"].iloc[i].split(";")
              res = list(filter(lambda x: "LOC" in x, t))
              if len(res) == 0:
                  res = list(filter(lambda x: "GeneID:" in x, t))
                  if len(res) != 0:
                      ID = res[0].split("GeneID:",1)[1]
                      ids.append(ID)
              else:
                  ID =res[0].split("LOC",1)[1]
                  ids.append(ID)
          #ids is a list of IDs we found in the rows where our loci is between start and stop
          uni=len(Counter(ids))

          if uni > 1:
              #if there is more than 1 kind of gene ID, there are overlapping genes, investigate cases separately
              new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "genes_overlap", "gene": "NA"},
                          ignore_index = True)

          gene_name = ids[0] #technically it might not be a gene

          if uni == 1:
              #if there is only one kind of region
              #checking what kind of region
              if temp['type'].str.contains('gene').any():
                if temp['type'].str.contains('pseudogene').any():
                      new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "pseudogene", "gene":"NA" },ignore_index = True)
                # snRNA, tRNA, lnc_RNA, snoRNA, rRNA, miRNA also have line with "gene"
                elif temp['type'].str.contains('snRNA').any():
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "snRNA", "gene":ids[0] },ignore_index = True)
                elif temp['type'].str.contains('tRNA').any():
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "tRNA", "gene":ids[0] },ignore_index = True)
                elif temp['type'].str.contains('lnc_RNA').any():
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "lnc_RNA", "gene":ids[0] },ignore_index = True)
                elif temp['type'].str.contains('snoRNA').any():
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "snoRNA", "gene":ids[0] },ignore_index = True)
                elif temp['type'].str.contains('rRNA').any():
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "rRNA", "gene":ids[0] },ignore_index = True)
                elif temp['type'].str.contains('miRNA').any():
                    new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "miRNA", "gene":ids[0] },ignore_index = True)
                else:
                  if temp['type'].str.contains('exon').any():
                      if temp['type'].str.contains('CDS').any():
                          #in exon
                          new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "exon", "gene":f"LOC{ids[0]}" },ignore_index = True)
                      else:
                          direction=temp['s2'].iloc[0]
                          e_start=temp[temp["type"]=="exon"].iloc[0].start
                          e_stop=temp[temp["type"]=="exon"].iloc[0].stop
                          g_start=temp[temp["type"]=="gene"].iloc[0].start
                          g_stop=temp[temp["type"]=="gene"].iloc[0].stop
                          if direction == "+":
                              if e_start == g_start:
                                  #print('5prime1')
                                  new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "5'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                              elif e_stop == g_stop:
                                  #print("3prime1")
                                  new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "3'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                              else:
                                  #print("ERROR1") #alternative UTR, inspect separately
                                  #print(temp)
                                  new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "alternative UTR", "gene":"NA" },ignore_index = True)
                          elif direction == "-":
                              if e_start == g_start:
                                  #print('3prime2')
                                  new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "3'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                              elif e_stop == g_stop:
                                  #print("5prime2")
                                  new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "5'UTR", "gene":f"LOC{ids[0]}" },ignore_index = True)
                              else:
                                  #print("ERROR2") #alternative UTR, inspect separately
                                  new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "alternative UTR", "gene":"NA" },ignore_index = True)
                          else:
                              print("ERROR3") #this should never return

                  else:
                      #in intron
                      new_df = new_df.append({'chr' : ch, 'pos' : pos, 'region' : "intron", "gene":f"LOC{ids[0]}" },ignore_index = True)

              else:
                  print("not gene")
  return new_df

def process_overlap(annotation_df,overl_df):

    new_df2=pd.DataFrame() #similar to code above

    for index, row in overl_df.iterrows():
        ch = row.chr
        pos = int(row.pos)
        temp=annotation_df[(annotation_df['chr']==ch) & (annotation_df['start']<=pos) & (annotation_df['stop']>=pos)]
        if len(temp) == 0:
            #print("not in annotation")
            new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "not_annot", "gene" : "NA"},
                            ignore_index = True)
        if len(temp) == 1:
            print("only one")
            print(temp['type'])
        if len(temp) > 1:
            #print("more than 1")
            #checking for gene overlap
            ids=[]
            for i in range(len(temp)):
                t=temp["info"].iloc[i].split(";")
                res = list(filter(lambda x: "LOC" in x, t))
                if len(res) == 0:
                    res = list(filter(lambda x: "GeneID:" in x, t))
                    if len(res) != 0:
                        ID = res[0].split("GeneID:",1)[1]
                        ids.append(ID)
                else:
                    ID =res[0].split("LOC",1)[1]
                    ids.append(ID)
            uni=len(Counter(ids))
            all_loc=Counter(ids)

            for i in all_loc:
                gt=annotation_df[(annotation_df['chr']==ch) & (annotation_df['start']<=pos) & (annotation_df['stop']>=pos) & (annotation_df['info'].str.contains(i)) ]
                if gt['type'].str.contains('gene').any():
                    if gt['type'].str.contains('pseudogene').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "pseudogene", "gene":"NA" },ignore_index = True)
                    # snRNA, tRNA, lnc_RNA, snoRNA, rRNA, miRNA also have line with "gene"
                    elif gt['type'].str.contains('snRNA').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "snRNA", "gene":i },ignore_index = True)
                    elif gt['type'].str.contains('tRNA').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "tRNA", "gene":i },ignore_index = True)
                    elif gt['type'].str.contains('lnc_RNA').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "lnc_RNA", "gene":i },ignore_index = True)
                    elif gt['type'].str.contains('snoRNA').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "snoRNA", "gene":i },ignore_index = True)
                    elif gt['type'].str.contains('rRNA').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "rRNA", "gene":i },ignore_index = True)
                    elif gt['type'].str.contains('miRNA').any():
                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "miRNA", "gene":i },ignore_index = True)
                    else:
                        if gt['type'].str.contains('exon').any():
                            if gt['type'].str.contains('CDS').any():
                                #print("in exon")
                                new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "exon", "gene":f"LOC{i}" },ignore_index = True)
                            else:
                                direction=gt['s2'].iloc[0]
                                e_start=gt[gt["type"]=="exon"].iloc[0].start
                                e_stop=gt[gt["type"]=="exon"].iloc[0].stop
                                g_start=gt[gt["type"]=="gene"].iloc[0].start
                                g_stop=gt[gt["type"]=="gene"].iloc[0].stop
                                if direction == "+":
                                    if e_start == g_start:
                                        #print('5prime1')
                                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "5'UTR", "gene":f"LOC{i}" },ignore_index = True)
                                    elif e_stop == g_stop:
                                        #print("3prime1")
                                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "3'UTR", "gene":f"LOC{i}" },ignore_index = True)
                                    else:
                                        #print("ERROR1") #alternative UTR?
                                        #print(gt)
                                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "alternative UTR", "gene":f"LOC{i}" },ignore_index = True)
                                elif direction == "-":
                                    if e_start == g_start:
                                        #print('3prime2')
                                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "3'UTR", "gene":f"LOC{i}" },ignore_index = True)
                                    elif e_stop == g_stop:
                                        #print("5prime2")
                                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "5'UTR", "gene":f"LOC{i}" },ignore_index = True)
                                    else:
                                        #print("ERROR2")
                                        new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "alternative UTR", "gene":f"LOC{i}" },ignore_index = True)
                                else:
                                    print("ERROR3")
                                    print(temp)
                                    #print(gt)

                        else:
                            #print("in intron")
                            new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "intron", "gene":f"LOC{i}" },ignore_index = True)

                else:
                    #print(temp)
                    #print("not gene")
                    new_df2 = new_df2.append({'chr' : ch, 'pos' : pos, 'region' : "pseudogene", "gene":"NA" },ignore_index = True)

    return new_df2

def process_annotation_data(annotation_df):

    promoters=annotation_df.copy()
    promoters.drop(promoters[promoters['type'] == "region"].index, inplace = True)

    imp_list=["gene","snRNA", "tRNA", "lnc_RNA", "snoRNA", "rRNA", "miRNA"] #regions in front of these regions will be labelled as promoter
    for l in imp_list:
        promoters.loc[(promoters["type"]==l) & (promoters["s2"]=="+"), "start"] -= 5000
        promoters.loc[(promoters["type"]==l) & (promoters["s2"]=="-"), "stop"] += 5000

    promoters['info'] = promoters['info'].str.replace(',',';')

    return promoters

def process_promoter(missing_annot,promoters):

    prom_df=pd.DataFrame()

    for index, row in missing_annot.iterrows(): #redo annotation but now annotation df is extended to include promoter regions (see process_annotation_data())
        ch = row.chr
        pos = int(row.pos)
        temp=promoters[(promoters['chr']==ch) & (promoters['start']<=pos) & (promoters['stop']>=pos)]
        if len(temp) == 0:
            #print("not in annotation")
            prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "not_annot", "gene" : "NA"},
                            ignore_index = True)
        else:
            ids=[]
            for i in range(len(temp)):
                t=temp["info"].iloc[i].split(";")
                res = list(filter(lambda x: "LOC" in x, t))
                if len(res) == 0:
                    res = list(filter(lambda x: "GeneID:" in x, t))
                    if len(res) != 0:
                        ID = res[0].split("GeneID:",1)[1]
                        #print(ID)
                        ids.append(ID)
                else:
                    ID =res[0].split("LOC",1)[1]
                    ids.append(ID)
            uni=len(Counter(ids))
            gene_name = ids[0]
            if uni > 1:
                #print("genes overlap")
                prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "genes_overlap", "gene": "NA"},
                            ignore_index = True)

            if uni == 1:
                #print("new")
                if temp['type'].str.contains('gene').any():
                    if temp['type'].str.contains('pseudogene').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "pseudogene", "gene":"NA" },ignore_index = True)
                    # snRNA, tRNA, lnc_RNA, snoRNA, rRNA, miRNA also have line with "gene"
                    elif temp['type'].str.contains('snRNA').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "snRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('tRNA').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "tRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('lnc_RNA').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "lnc_RNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('snoRNA').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "snoRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('rRNA').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "rRNA", "gene":ids[0] },ignore_index = True)
                    elif temp['type'].str.contains('miRNA').any():
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "miRNA", "gene":ids[0] },ignore_index = True)
                    else:
                        prom_df = prom_df.append({'chr' : ch, 'pos' : pos, 'region' : "protein_coding", "gene" : f"LOC{gene_name}"},
                                        ignore_index = True)

    return prom_df

def check_enhancers(only_not_annot_prom, atac, chip, lvar,lnc): # NOTE: didn't end up using it
    #takes loci not annotated, outputs df of loci that fall in ATAC-seq, Chip-seq, etc. regions

    enh_check=pd.DataFrame()

    for index, row in only_not_annot_prom.iterrows():
        ch = row.chr
        pos = row.pos

        temp=atac[(atac['chr']==ch) & (atac['start']<=pos) & (atac['stop']>=pos)]
        if len(temp) > 0:
            #print("in atac")
            enh_check = enh_check.append({'chr' : ch, 'pos' : pos, 'region' : "atac" },ignore_index = True)

        temp=chip[(chip['chr']==ch) & (chip['start']<=pos) & (chip['stop']>=pos)]
        if len(temp) > 0:
            #print("in chip")
            enh_check = enh_check.append({'chr' : ch, 'pos' : pos, 'region' : "chip" },ignore_index = True)

        temp=lvar[(lvar['chr']==ch) & (lvar['start']<=pos) & (lvar['stop']>=pos)]
        if len(temp) > 0:
            #print("in lvar")
            enh_check = enh_check.append({'chr' : ch, 'pos' : pos, 'region' : "lvar" },ignore_index = True)

        temp=lnc[(lnc['chr']==ch) & (lnc['start']<=pos) & (lnc['stop']>=pos)]
        if len(temp) > 0:
            #print("in lnc")
            enh_check = enh_check.append({'chr' : ch, 'pos' : pos, 'region' : "lnc" },ignore_index = True)

    return enh_check
