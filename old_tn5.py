import Levenshtein




def judge_adpter(str1):
    f1="ATGATGTGTTTGAGTATTGTTTGAGATGTGTATAAGAGATAG"

    my_seq=str1.replace("C","T")


    for i in range(len(my_seq)-15 +1):
        a1=my_seq[i:i+15]
        for j in range(len(f1)-15 +1):
            a2=f1[j:j+15]
            if Levenshtein.distance(a1,a2)<=2:
                return 1

    return 0







c1=[]
c2=[]

with open("cgs_CGS-Tn-ES_1.fq")as file1,open("cgs_CGS-Tn-ES_2.fq")as file2,\
        open("garbage_1.fq","w")as file3,open("garbage_2.fq","w")as file4,\
        open("noadpter_1.fq","w")as file5,open("noadpter_2.fq","w")as file6,\
        open("ATGATGTG_1.fq","w")as file7,open("ATGATGTG_2.fq","w")as file8:


    for i,line in enumerate(zip(file1,file2)):


        if i%10000==0:
            print(i)

        if i%4==0:
            name1=line[0].strip()
            name2=line[1].strip()
        if i %4==1:
            read1=line[0].strip()
            read2=line[1].strip()

            if read1[-10:]=="G"*10 or read2[-10:]=="G"*10:
                print(name1,file=file3)
                print(read1,file=file3)
                print("+",file=file3)
                print("I"*len(read1),file=file3)

                print(name2,file=file4)
                print(read2,file=file4)
                print("+",file=file4)
                print("I"*len(read2),file=file4)


            else:
                f_score=judge_adpter(read1[:35])
                if f_score==1:
                    print(name1,file=file7)
                    print(read1,file=file7)
                    print("+",file=file7)
                    print("I"*len(read1),file=file7)

                    print(name2,file=file8)
                    print(read2,file=file8)
                    print("+",file=file8)
                    print("I"*len(read2),file=file8)

                else:
                    print(name1, file=file5)
                    print(read1, file=file5)
                    print("+", file=file5)
                    print("I" * len(read1), file=file5)

                    print(name2, file=file6)
                    print(read2, file=file6)
                    print("+", file=file6)
                    print("I" * len(read2), file=file6)







