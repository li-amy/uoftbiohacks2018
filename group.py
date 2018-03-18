from setUpCh import chSetUp

# write gene.json file

def geneJsonWrite(chDF, fName):
    funcIntx = dict()
    with open("Chr20FuncIntx.tsv") as i:
        first_line = i.readline()
        print(first_line)
        for line in i:
            values = line.strip().split('\t')
            if values[0] not in funcIntx:
                funcIntx[values[0]] = [(values[1], values[2])]
            else:
                funcIntx[values[0]].append((values[1], values[2]))
    # print(funcIntx)
    with open(fName, 'w') as o:
        o.write("[\n")
        for g in chDF.index:
            gene_type = chDF.loc[g].type
            goa_d = chDF.loc[g].GOAdescr
            imports = ""
            if g not in funcIntx:
                imports = ""
            else:
                if len(funcIntx[g]) == 1:
                    ig = funcIntx[g][0][0]
                    itype = chDF.loc[ig].type
                    igoad = chDF.loc[ig].GOAdescr
                    imports += "\"gene.%s.%s.%s\"" % (itype,igoad, ig)
                else:
                    for i in range(len(funcIntx[g])):
                        ig = funcIntx[g][i][0]
                        itype = chDF.loc[ig].type
                        igoad = chDF.loc[ig].GOAdescr
                        imports += "\"gene.%s.%s.%s\"," % (itype,igoad, ig)
                    imports = imports[:-1]
            # imports = "gene.protein_coding.ribonucleoprotein complex assembly.AAR2"
            line = "{\"name\":\"gene.%s.%s.%s\",\"size\":0,\"imports\":[%s]},\n" % (gene_type, goa_d, g, imports)
            o.write(line)
        o.write("]\n")

    unique_goa_descr = dict()
    unique_gene_type = dict()
    for g in chDF.index:
        # if g == "RALGAPB":
        #     print(chDF.loc[g])
        if chDF.loc[g].GOAdescr not in unique_goa_descr:
            unique_goa_descr[chDF.loc[g].GOAdescr] = 1
        else:
            unique_goa_descr[chDF.loc[g].GOAdescr] += 1
        if chDF.loc[g].type not in unique_gene_type:
            unique_gene_type[chDF.loc[g].type] = 1
        else:
            unique_gene_type[chDF.loc[g].type] += 1
    # print(unique_goa_descr["ribonucleoprotein complex assembly"])
    # print(unique_gene_type)
    # print(len(unique_gene_type))

if __name__ == '__main__':
    print("setting up...")
    ch20 = chSetUp()
    geneJsonWrite(ch20, "gene.json")
