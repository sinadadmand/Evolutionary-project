from Bio import SeqIO
from math import log2
from numpy.core import float32
import time


def fasta_to_entry_pit(infile, outfile):
    starttime = time.time()  # timer start
    header = "Entry,Organism_ID,Pit,len,A,C,D,E,F,G,H,I,K,L,M,N,O,P,Q,R,S,T,U,V,W,Y\n"  # output file header
    with open(outfile, "w") as filewrite:
        filewrite.write(header)  # write header
        for record in SeqIO.parse(infile, "fasta"):  # input Fasta file
            firstchEntry = record.id.index("|", 1)  # export Entry from Fasta
            lastchEntry = record.id.index("|", 5)  # export Entry from Fasta
            firstchOrg = record.description.index(" OX=", -170)  # export Organism ID from Fasta
            lastchOrg = record.description.index(" ", firstchOrg + 5)  # export Organism ID from Fasta
            Organism_ID = record.description[firstchOrg + 4:lastchOrg]  # write Organism to file
            Entry = record.id[firstchEntry + 1:lastchEntry]  # write Entry to file
            Sequence = record.seq  # write Sequence to file

            A = (Sequence.count("A"))  # count A string in Sequence

            C = (Sequence.count("C"))  # count C string in Sequence
            D = (Sequence.count("D") + (Sequence.count("B") / 2))  # count D string in Sequence
            E = (Sequence.count("E") + (Sequence.count("Z") / 2))  # count E string in Sequence
            F = (Sequence.count("F"))  # count F string in Sequence
            G = (Sequence.count("G"))  # count G string in Sequence
            H = (Sequence.count("H"))  # count H string in Sequence
            I = (Sequence.count("I") + (Sequence.count("J") / 2))  # count I string in Sequence

            K = (Sequence.count("K"))  # count K string in Sequence
            L = (Sequence.count("L") + (Sequence.count("J") / 2))  # count L string in Sequence
            M = (Sequence.count("M"))  # count M string in Sequence
            N = (Sequence.count("N") + (Sequence.count("B") / 2))  # count N string in Sequence
            O = (Sequence.count("O"))  # count O string in Sequence
            P = (Sequence.count("P"))  # count P string in Sequence
            Q = (Sequence.count("Q") + (Sequence.count("Z") / 2))  # count Q string in Sequence
            R = (Sequence.count("R"))  # count R string in Sequence
            S = (Sequence.count("S"))  # count S string in Sequence
            T = (Sequence.count("T"))  # count T string in Sequence
            U = (Sequence.count("U"))  # count U string in Sequence
            V = (Sequence.count("V"))  # count V string in Sequence
            W = (Sequence.count("W"))  # count W string in Sequence

            Y = (Sequence.count("Y"))  # count Y string in Sequence

            length = int(A + C + D + E + F + G + H + I + K + L + M + N + O + P + Q + R + S + T + U + V + W + Y)  #
            # sequence length

            Pit = float32(  # calculate Pit
                log2((A / length) ** (-A / length)) + log2((C / length) ** (-C / length)) + log2(
                    (D / length) ** (-D / length)) + log2(
                    (E / length) ** (-E / length)) + log2((F / length) ** (-F / length)) + log2(
                    (G / length) ** (-G / length)) + log2(
                    (H / length) ** (-H / length)) + log2((I / length) ** (-I / length)) + log2(
                    (K / length) ** (-K / length)) + log2(
                    (L / length) ** (-L / length)) + log2((M / length) ** (-M / length)) + log2(
                    (N / length) ** (-N / length)) + log2(
                    (O / length) ** (-O / length)) + log2((P / length) ** (-P / length)) + log2(
                    (Q / length) ** (-Q / length)) + log2(
                    (R / length) ** (-R / length)) + log2((S / length) ** (-S / length)) + log2(
                    (T / length) ** (-T / length)) + log2(
                    (U / length) ** (-U / length)) + log2((V / length) ** (-V / length)) + log2(
                    (W / length) ** (-W / length)) + log2((Y / length) ** (-Y / length)))

            output = str(Entry) + "," + str(Organism_ID) + "," + str(Pit) + "," + str(length) + "," + str((
                A)) + "," + str((C)) + "," + str((D)) + "," + str((E)) + "," + str((F)) + "," + str(
                (G)) + "," + str((H)) + "," + str((
                I)) + "," + str((K)) + "," + str((L)) + "," + str((M)) + "," + str((N)) + "," + str(
                (O)) + "," + str((P)) + "," + str((
                Q)) + "," + str((R)) + "," + str((S)) + "," + str((T)) + "," + str((U)) + "," + str(
                (V)) + "," + str((W)) + "," + str((
                Y)) + "\n"
            filewrite.write(output)
    endtime = time.time()
    timex = endtime - starttime
    print("total time taken is {:3.0f} hours {:3.0f} minutes {:3.0f} seconds {:3.0f} milisec".format((timex // 3600), (
                timex % 3600) // 60, (timex % 60 // 1), (timex - int(timex))))
