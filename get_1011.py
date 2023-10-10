import Levenshtein
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Align.substitution_matrices import Array
from Bio.Align import PairwiseAligner
def read_fq(fq_file):
    fq_dir={}
    with open(fq_file)as file:
        for i,line in enumerate(file):
            if i%4==0:
                name=line.split()[0]
            if i%4==1:
                seq=line.strip()
            if i%4==3:
                q=line.strip()
                fq_dir[name[1:]]=(seq,q)
    return fq_dir


fq1=read_fq("/home/wangzc/hair/data/trimmed/CS-P0_80_90_1.fq")
fq2=read_fq("/home/wangzc/hair/data/trimmed/CS-P0_80_90_2.fq")


def readgenome(name):
    dd = {}
    with open(name, 'r') as input_fasta:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            dd[record.id] = str(record.seq)
    return dd

dd=readgenome("/data/wangzc/cgs/new_mm9_lam_mc.fa")


def complement_dna(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    reverse = dna[::-1].translate(complement)
    return reverse
def check(str1, str2):
    length1 = len(str1)
    length = min(length1, len(str2))
    k = max(range(length + 1),
            key=lambda i: i if Levenshtein.hamming(str1[length1 - i:], str2[:i]) < i * 0.1 else False)
    return k

def trim_overlap(read1, read2):
    a1 = check(read1, complement_dna(read2))
    if a1 > 5:
        return len(read1)+len(read2)-a1
    else:
        a2 = check(complement_dna(read2), read1)
        if a2 > min(len(read2), len(read1)) * 0.8 and a2 > 5:
            return 900

    return  1000



class Alignment:
    alphabet = "ACGTN"
    rule_matrix_1 = np.array([
        [1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    rule_matrix_2 = np.array([
        [1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    def __init__(self, sub_matrix, mode='global', open_gap_score=-5, extend_gap_score=-2):
        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score
        if sub_matrix == 'rule_matrix_1':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_1)
        elif sub_matrix == 'rule_matrix_2':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_2)
        else:
            raise ValueError(f"Invalid sub_matrix argument: {sub_matrix}")
        self.aligner = aligner
        self.type = sub_matrix

    def _complement_dna(self,dna):
        complement = str.maketrans('ATCG', 'TAGC')
        reverse = dna[::-1].translate(complement)
        return reverse

    def align(self, seq1, seq2):
        if len(seq1)==0 or len(seq2)==0:
            return "N", "N", 1



        alignments = self.aligner.align(self._complement_dna(seq1),self._complement_dna(seq2))
        align_read1, align_read2 = self._complement_dna(alignments[0][0]), self._complement_dna(alignments[0][1])
        if self.type == 'rule_matrix_1':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): -10, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
            }

            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1




        elif self.type == 'rule_matrix_2':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): 1, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10,  ('-', '-'): -10
            }
            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1


        else:
            raise ValueError(f"Invalid sub_matrix argument: {self.type}")

        return  align_read1, align_read2,max_index




align_2=Alignment('rule_matrix_2')
align_1=Alignment('rule_matrix_1')

def jug_mode2(li_left,li_right,li_left_yz,li_right_yz):   # 这个函数保证了结果，max_length>=3,max_start_index<=3
    re_li_right=complement_dna(li_right)
    yz1_len=0
    for i in range(min(len(li_left),len(li_right))):
        if li_left[i]==re_li_right[i]:
            yz1_len=yz1_len+1
        else:
            break

    re_li_left_yz=complement_dna(li_left_yz)
    re_yz_len=0
    for i in range(min(len(re_li_left_yz),len(li_right_yz))):
        if re_li_left_yz[i]==li_right_yz[i]:
            re_yz_len=re_yz_len+1
        else:
            break


    if re_yz_len == 0 and yz1_len>=3:
        max_length=yz1_len
        max_start_index=0

    elif re_yz_len!=0 and yz1_len+re_yz_len>=3:
        max_length = yz1_len + re_yz_len
        max_start_index = -re_yz_len
    else:                    #(re_yz_len!=0 and yz1_len+re_yz_len<3 ) or (re_yz_len==0 and yz1_len<3)
        # 这里需要执行，磁铁起点向后位置1-3位
        max_length = 0  # 最长匹配子串的长度
        max_start_index = 0  # 最长匹配子串的起始索引
        current_length = 0  # 当前匹配子串的长度
        current_start_index = 0  # 当前匹配子串的起始索引
        for i in range(min(len(li_left),len(li_right))):
            if li_left[i] == re_li_right[i]:
                current_length += 1
                if current_length > max_length:
                    max_length = current_length
                    max_start_index = current_start_index
            else:
                current_length = 0
                current_start_index = i + 1
        if max_start_index == 0:
            re_yz_str1 = complement_dna(li_left_yz)
            re_len = 0
            for i in range(min(len(li_left_yz),len(li_right_yz))):
                if re_yz_str1[i] == li_right_yz[i]:
                    re_len = re_len + 1
                else:
                    break
            if re_len != 0:
                max_length = max_length + re_len
                max_start_index = -re_len
        # 执行完毕
        if max_length<3 or max_start_index>3:
            max_length=0
            max_start_index=0
    return max_length,max_start_index










def yz_mode3_windows(target_seq,ref_seq):
    my_list=[]
    target_len=len(target_seq)
    ref_start=-1
    for i in range(len(ref_seq)-target_len):
        temp_ref=ref_seq[i:i+target_len]
        temp_sc=Levenshtein.distance(target_seq.replace("C","T"),temp_ref.replace("C","T"))
        if temp_sc<=target_len*0.1:
            ref_start=i
            my_list.append((i,temp_sc))

    if len(my_list)!=0:
        result = min(my_list, key=lambda x: x[1])
        ref_start=result[0]


    return ref_start


def find_longest_matching_substring_mode3(str1, str2):
    max_length = 0  # 最长匹配子串的长度
    max_start_index = 0  # 最长匹配子串的起始索引
    current_length = 0  # 当前匹配子串的长度
    current_start_index = 0  # 当前匹配子串的起始索引

    for i in range(min(len(str1),len(str2))):
        if str1[i] == str2[i]:
            current_length += 1
            if current_length > max_length:
                max_length = current_length
                max_start_index = current_start_index
        else:
            current_length = 0
            current_start_index = i + 1

    return max_length,max_start_index
def jug_mode3(li_left,li_right,li_left_yz,li_right_yz):
    re_li_right=complement_dna(li_right)
    yz1_len=0
    for i in range(min(len(li_left),len(li_right))):
        if li_left[i]==re_li_right[i]:
            yz1_len=yz1_len+1
        else:
            break

    re_li_left_yz=complement_dna(li_left_yz)
    re_yz_len=0
    for i in range(min(len(re_li_left_yz),len(li_right_yz))):
        if re_li_left_yz[i]==li_right_yz[i]:
            re_yz_len=re_yz_len+1
        else:
            break

    if re_yz_len == 0 and yz1_len>=3:
        max_length=yz1_len
        max_start_index=0

    elif re_yz_len!=0 and yz1_len+re_yz_len>=3:
        max_length = yz1_len + re_yz_len
        max_start_index = -re_yz_len
    else:                    #(re_yz_len!=0 and yz1_len+re_yz_len<3 ) or (re_yz_len==0 and yz1_len<3)

        max_length, max_start_index = find_longest_matching_substring_mode3(li_right_yz[::],complement_dna(li_left_yz[::]) )


        if max_length<3 or max_start_index>3:
            max_length=0
            max_start_index=0
        else:
            max_length=max_length
            max_start_index=-(max_length+max_start_index)

    return max_length,max_start_index









def my_fun(fq1_seq,fq2_seq,ref_seq):
    seq1=fq1_seq[0]
    seq2=fq2_seq[0]
    q1=fq1_seq[1]
    q2=fq2_seq[1]

    # seq1=fq1_seq
    # seq2=fq2_seq



    a_seq1,a_seq2,len1 = align_2.align(seq1, seq2)
    mode2_index=0
    mode3_index=0
    mode4_index=0




    # 执行模式2
    jug_len = trim_overlap(seq1, seq2) # 这里给出更严格overlap限制，最短为11
    if jug_len < 290:
        mode2=1
        if jug_len % 2 == 0:
            half_seq2 = complement_dna(seq2[:int(jug_len / 2)])
            half_ref = ref_seq[int(jug_len / 2) :int(jug_len / 2)  + int(jug_len / 2)]
            h_seq1, h_seq2, len2 = align_1.align(half_seq2, half_ref)
            li_right = ref_seq[int(jug_len / 2) :int(jug_len / 2)  + len2]
            li_left = ref_seq[int(jug_len / 2) - len2:int(jug_len / 2)]
            li_left_yz = ref_seq[int(jug_len / 2) - len2-5:int(jug_len / 2) - len2]
            li_right_yz = ref_seq[int(jug_len / 2) + len2:int(jug_len / 2) + len2+5]
            mode2_len,mode2_index=jug_mode2(li_left,li_right,li_left_yz,li_right_yz)


            li_left_index =int(jug_len / 2) - len2
            li_right_index=int(jug_len / 2)  + len2

        else:
            half_seq2 = complement_dna(seq2[:int(jug_len / 2)])
            half_ref = ref_seq[int(jug_len / 2) + 1:int(jug_len / 2) + 1 + int(jug_len / 2)]
            h_seq1, h_seq2, len2 = align_1.align(half_seq2, half_ref)
            li_right=ref_seq[int(jug_len / 2) + 1:int(jug_len / 2) + 1+len2]
            li_left = ref_seq[int(jug_len / 2) - len2:int(jug_len / 2) ]
            li_left_yz=ref_seq[int(jug_len / 2) - len2-5:int(jug_len / 2) - len2 ]
            li_right_yz = ref_seq[int(jug_len / 2) + 1 + len2 :int(jug_len / 2) + 1 + len2 +5 ]
            mode2_len,mode2_index=jug_mode2(li_left,li_right,li_left_yz,li_right_yz)

            li_left_index =int(jug_len / 2) - len2
            li_right_index=int(jug_len / 2) + 1+len2


        if mode2_len >= 3:
            ci_huan_mode2 = ref_seq[li_left_index + mode2_index:li_right_index - mode2_index]
            ci_wai_yz_mode2 = ref_seq[li_left_index + mode2_index - 10:li_left_index + mode2_index] + " " + \
                        ref_seq[ li_right_index - mode2_index:li_right_index - mode2_index + 10]
            ci_len_mode2 = mode2_len
        else:
            ci_huan_mode2="N"
            ci_wai_yz_mode2="N"
            ci_len_mode2=0


    elif jug_len==900:
        mode2 = 2
        jug_len = len(seq1)+len(seq2)-check(complement_dna(seq2), seq1)

        if jug_len % 2 == 0:
            temp_len_9=len(seq2)-int(jug_len / 2)
            half_seq2 = complement_dna(seq2[:temp_len_9])
            half_ref = ref_seq[temp_len_9 : temp_len_9 + temp_len_9]
            h_seq1, h_seq2, len2 = align_1.align(half_seq2, half_ref)
            li_right = ref_seq[temp_len_9 :temp_len_9  + len2]
            li_left = ref_seq[temp_len_9 - len2:temp_len_9]
            li_left_yz = ref_seq[temp_len_9 - len2-5:temp_len_9 - len2]
            li_right_yz = ref_seq[temp_len_9 + len2:temp_len_9 + len2+5]

            mode2_len,mode2_index=jug_mode2(li_left,li_right,li_left_yz,li_right_yz)

            li_left_index =temp_len_9 - len2
            li_right_index=temp_len_9  + len2
        else:
            temp_len_9 = len(seq2) - int(jug_len / 2) -1
            half_seq2 = complement_dna(seq2[:temp_len_9])
            half_ref = ref_seq[temp_len_9 + 1:temp_len_9 + 1 + temp_len_9]
            h_seq1, h_seq2, len2 = align_1.align(half_seq2, half_ref)
            li_right=ref_seq[temp_len_9 + 1:temp_len_9 + 1+len2]
            li_left = ref_seq[temp_len_9 - len2:temp_len_9 ]
            li_left_yz=ref_seq[temp_len_9 - len2-5:temp_len_9 - len2 ]
            li_right_yz = ref_seq[temp_len_9 + 1 + len2 :temp_len_9 + 1 + len2 +5 ]

            mode2_len,mode2_index=jug_mode2(li_left,li_right,li_left_yz,li_right_yz)

            li_left_index =temp_len_9 - len2
            li_right_index=temp_len_9 + 1+len2


        if mode2_len >= 3:
            ci_huan_mode2 = ref_seq[li_left_index + mode2_index:li_right_index - mode2_index]
            ci_wai_yz_mode2 = ref_seq[li_left_index + mode2_index - 10:li_left_index + mode2_index] + " " + \
                        ref_seq[ li_right_index - mode2_index:li_right_index - mode2_index + 10]
            ci_len_mode2 = mode2_len
        else:
            ci_huan_mode2="N"
            ci_wai_yz_mode2="N"
            ci_len_mode2=0

    else:
        mode2 = 3
        ci_huan_mode2 = "N"
        ci_wai_yz_mode2 = "N"
        ci_len_mode2 = 0








    # 执行mode1，这里mode1_len，是环磁的右终点，   len1为左茎+左磁铁，   磁铁是mode1_len-150+mode1_len_2
    # （这里mode1_len_2是右茎一部分+右磁铁），所以mode1_len+mode1_len_2为总长度多一个磁铁


    mode1_seq1, mode1_seq2, mode1_len_1 = align_1.align(seq1, ref_seq[:150])
    if mode1_len_1<145 and len1<140:
        mode1=1
        mode1_seq1_2, mode1_seq2_2, mode1_len_2 = align_2.align(seq1[len1:][::-1], seq2[len1:][::-1])

        if len1 - (mode1_len_1- 150+ mode1_len_2)<=0 or len1 - (mode1_len_1- 150+ mode1_len_2) >=mode1_len_1:
            mode1 = 0
            ci_huan_mode1 = "N"
            ci_wai_yz_mode1 = "N"
            ci_len_mode1 = 0


        else:
            li_left=ref_seq[len1 - (mode1_len_1- 150+ mode1_len_2): mode1_len_1]
            li_right=ref_seq[len1 - (mode1_len_1- 150+ mode1_len_2): mode1_len_1]
            li_left_yz=ref_seq[len1 - (mode1_len_1- 150+ mode1_len_2)-5: len1 - (mode1_len_1- 150+ mode1_len_2)]
            li_right_yz=ref_seq[mode1_len_1: mode1_len_1+5]

            mode1_len, mode1_index = jug_mode2(li_left, li_right, li_left_yz, li_right_yz)

            li_left_index =len1 - (mode1_len_1- 150+ mode1_len_2)
            li_right_index=mode1_len_1


            if mode1_len >= 3:
                ci_huan_mode1 = ref_seq[li_left_index + mode1_index:li_right_index - mode1_index]
                ci_wai_yz_mode1 = ref_seq[li_left_index + mode1_index - 10:li_left_index + mode1_index] + " " + \
                                  ref_seq[li_right_index - mode1_index:li_right_index - mode1_index + 10]
                ci_len_mode1 = mode1_len
            else:
                ci_huan_mode1 = "N"
                ci_wai_yz_mode1 = "N"
                ci_len_mode1 = 0

    else:
        mode1=0
        ci_huan_mode1 = "N"
        ci_wai_yz_mode1 = "N"
        ci_len_mode1 = 0





    # 开始执行模式3 这里需要限制用于比对的最小长度，那就10/11吧
    if 150-len1>10:
        mode3=1

        mode3_target=complement_dna(seq2[len1:])
        mode3_ref=ref_seq[len1:350]
        yz_start=yz_mode3_windows(mode3_target,mode3_ref)
        if yz_start< 0:
            ci_huan_mode3 = "N"
            ci_wai_yz_mode3 = "N"
            ci_len_mode3 = 0
        else:
            temp_huan=ref_seq[len1:yz_start+len(seq2)]
            yz_mode3_left=ref_seq[len1-20:len1]
            yz_mode3_right=ref_seq[yz_start + len(seq2):yz_start + len(seq2)+20]
            mode3_len, mode3_index = jug_mode3(temp_huan, temp_huan, yz_mode3_left, yz_mode3_right)


            if mode3_len<3:
                ci_huan_mode3 = "N"
                ci_wai_yz_mode3 = "N"
                ci_len_mode3 = 0
            else:
                ci_huan_mode3 = ref_seq[len1+mode3_index:yz_start+len(seq2)-mode3_index]
                ci_wai_yz_mode3=ref_seq[len1+mode3_index-10:len1+mode3_index]+" "+ref_seq[yz_start+len(seq2)-mode3_index:yz_start+len(seq2)-mode3_index+10]
                ci_len_mode3 = mode3_len

    else:
        mode3=0
        ci_huan_mode3 = "N"
        ci_wai_yz_mode3 = "N"
        ci_len_mode3 = 0



    # 开始执行模式4,对模式3的严谨补充，这里需要限制用于比对的最小长度，那就10/11吧
    if 150 - len1 > 10:
        mode4 = 1
        mode4_remove_len=0
        for i in range(140):
            if seq1[-i]==seq2[-i] or (seq1[-i]=="T" and seq2[-i]=="C") or (seq1[-i]=="G" and seq2[-i]=="A"):
                mode4_remove_len=mode4_remove_len+1
            else:
                break
        prime_len_mode4=int((len(seq2)-len1)*0.2)
        if mode4_remove_len<prime_len_mode4:
            mode4_remove_len=prime_len_mode4
        else:
            if (len(seq2)-len1)-mode4_remove_len>15:
                mode4_remove_len=mode4_remove_len-4
        mode4_target = complement_dna(seq2[len1:-mode4_remove_len])
        mode4_ref = ref_seq[len1:350]
        yz_start = yz_mode3_windows(mode4_target, mode4_ref)



        if yz_start < 0:
            ci_huan_mode4 = "N"
            ci_wai_yz_mode4 = "N"
            ci_len_mode4 = 0
        else:
            temp_huan = ref_seq[len1: yz_start + len(seq2) - mode4_remove_len]
            yz_mode4_left = ref_seq[len1 - 20: len1]
            yz_mode4_right = ref_seq[ yz_start + len(seq2) - mode4_remove_len: yz_start + len(seq2) - mode4_remove_len + 20]
            mode4_len, mode4_index = jug_mode3(temp_huan, temp_huan, yz_mode4_left, yz_mode4_right)


            if mode4_len<3:
                ci_huan_mode4 = "N"
                ci_wai_yz_mode4 = "N"
                ci_len_mode4 = 0
            else:
                ci_huan_mode4 = ref_seq[len1+mode4_index:yz_start+len(seq2)-mode4_remove_len-mode4_index]
                ci_wai_yz_mode4=ref_seq[len1+mode4_index-10:len1+mode4_index]+" "+ref_seq[yz_start+len(seq2)-mode4_remove_len-mode4_index:yz_start+len(seq2)-mode4_remove_len-mode4_index+10]
                ci_len_mode4 = mode4_len
    else:
        mode4 = 0
        ci_huan_mode4 = "N"
        ci_wai_yz_mode4 = "N"
        ci_len_mode4 = 0



    final_len_list=[ci_len_mode1,ci_len_mode2,ci_len_mode3,ci_len_mode4]
    max_index = final_len_list.index(max(final_len_list))
    if max_index==0:
        return ci_huan_mode1,ci_wai_yz_mode1,ci_len_mode1
    elif max_index==1:
        return ci_huan_mode2, ci_wai_yz_mode2, ci_len_mode2
    elif max_index==2:
        return ci_huan_mode3, ci_wai_yz_mode3, ci_len_mode3
    else:
        return ci_huan_mode4, ci_wai_yz_mode4, ci_len_mode4











pp=0


w_file=open("/home/wangzc/hair/data/trimmed/wwwwwwwwwwwwwwwwwwwwww.f","w")
samfile = pysam.AlignmentFile("/home/wangzc/hair/data/trimmed/unique_alignments.sam", "r")
for alignment in samfile:
    pp=pp+1
    print(pp)

    if pp%1000000==0:
        break
    if alignment.is_unmapped:
        continue


    ref_chr=alignment.reference_name
    read_name=alignment.qname

    if alignment.is_reverse:
        fq1_seq=fq1[read_name]
        fq2_seq=fq2[read_name]
        if alignment.reference_start+alignment.query_alignment_end-500<0:
            ref_seq = dd[ref_chr][:alignment.reference_start + alignment.query_alignment_end]
        else:
            ref_seq=dd[ref_chr][alignment.reference_start+alignment.query_alignment_end-500:alignment.reference_start+alignment.query_alignment_end]
        ref_seq=complement_dna(ref_seq)

        f_a,f_b,f_c=my_fun(fq1_seq, fq2_seq, ref_seq)
        print("@"+read_name,file=w_file)
        print(f_a,file=w_file)
        print(f_c,file=w_file)
        print(f_b,file=w_file)



    else:
        fq1_seq=fq1[read_name]
        fq2_seq=fq2[read_name]
        ref_seq = dd[ref_chr][alignment.reference_start:alignment.reference_start + 500]

        f_a,f_b,f_c=my_fun(fq1_seq, fq2_seq, ref_seq)
        print("@"+read_name,file=w_file)
        print(f_a,file=w_file)
        print(f_c,file=w_file)
        print(f_b,file=w_file)





samfile.close()






