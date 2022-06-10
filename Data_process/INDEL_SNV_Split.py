path_in = "/h/xutong/project_RV/21_nonSigGWAS/01_data/filter_V2/"
path_out = "/h/xutong/project_RV/21_nonSigGWAS/08_INDEL/01_preData/"

filenm = "04_HotSpots_ConsFilter.tsv"
file_in = path_in + filenm
file_INDEL = path_out + "TestData_INDEL.tsv"
file_SNV = path_out + "TestData_SNV.tsv"

FO_INDEL = open(file_INDEL, "a+")
FO_SNV = open(file_SNV, "a+")

with open(file_in, "r") as FI:
    lines = FI.readlines()[1:]
    for line in lines:
        line = line.strip("\n")
        line = line.split("\t")
        
        if (len(line[2]) > 50) or (len(line[3]) > 50):
            pass
        
        elif (len(line[2]) > 1) or (len(line[3]) > 1):
            ## missing
            if len(line[2]) > len(line[3]):
                
                begin_pos = int(line[1]) - 1
                end_pos = int(line[1]) + (abs(len(line[2]) - len(line[3])))
                pos_id = "_".join(line[0:4])
                
                for tempos in range(begin_pos, end_pos):
                    btemp_pos = tempos
                    etemp_pos = tempos + 1
                    out_line = line[0:1] + [str(btemp_pos), str(etemp_pos)] + [pos_id]
                    out_line = "\t".join(out_line)
                    
                    FO_INDEL.write(out_line +  "\n")
            ## insert      
            elif len(line[2]) < len(line[3]):
                begin_pos = int(line[1]) - 1
                end_pos = int(line[1]) + (abs(len(line[2]) - len(line[3])))
                pos_id = "_".join(line[0:4])
                
                for tempos in [begin_pos - 1, end_pos]:
                    btemp_pos = tempos
                    etemp_pos = tempos + 1
                    out_line = line[0:1] + [str(btemp_pos), str(etemp_pos)] + [pos_id]
                    out_line = "\t".join(out_line)
                    
                    FO_INDEL.write(out_line +  "\n")
            ## Other INDEL    
            else:
                pass
        
        ## SNVs
        else:
            pos_id = "_".join(line[0:4])
            out_line = line[0:4] + [pos_id]
            out_line = "\t".join(out_line)
            
            FO_SNV.write(out_line +  "\n")
        
        
FO_INDEL.close()
FO_SNV.close()
