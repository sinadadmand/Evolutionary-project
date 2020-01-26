import pandas as pd
from numpy import log2
import time


def entry_to_species(sprot_input_file='Entries sprot.csv', trembl_input_file='Entries trembl.csv',
                     outfile='Species.csv', chunksize=8000000, merge='True'):
    starttime = time.time()
    out = pd.DataFrame()

    data = pd.read_csv(trembl_input_file, chunksize=chunksize)

    for chunk in data:
        chunk = chunk.drop(columns=['Pit'])
        group = chunk.groupby('Organism_ID').sum()
        out = out.append(group)

    if merge == 'True':
        data2 = pd.read_csv(sprot_input_file).drop(columns=['Pit'])
        data2 = data2.groupby('Organism_ID').sum()
        out = pd.concat([out, data2])

    out = out.groupby('Organism_ID').sum()
    out.to_csv(outfile)
    endtime = time.time()
    timex = endtime - starttime
    print(
        "total time taken is {:3.0f} hours {:3.0f} minutes {:3.0f} seconds {:3.0f} milisecs".format((timex // 3600), (
                timex % 3600) // 60,
                                                                                                    (timex % 60 // 1),
                                                                                                    (timex - int(
                                                                                                        timex))))


def species_pit(input_file='Species.csv', outfile='Species Pit.csv', proteomes_file="proteomes-redundant-no.tab"):
    starttime = time.time()
    data = pd.DataFrame(pd.read_csv(input_file, dtype=float))
    length = data.len
    Pit = pd.DataFrame(log2((data["A"] / length) ** (-data["A"] / length)) + log2(
        (data["C"] / length) ** (-data["C"] / length)) + log2(
        (data["D"] / length) ** (-data["D"] / length)) + log2(
        (data["E"] / length) ** (-data["E"] / length)) + log2((data["F"] / length) ** (-data["F"] / length)) + log2(
        (data["G"] / length) ** (-data["G"] / length)) + log2(
        (data["H"] / length) ** (-data["H"] / length)) + log2((data["I"] / length) ** (-data["I"] / length)) + log2(
        (data["K"] / length) ** (-data["K"] / length)) + log2(
        (data["L"] / length) ** (-data["L"] / length)) + log2((data["M"] / length) ** (-data["M"] / length)) + log2(
        (data["N"] / length) ** (-data["N"] / length)) + log2(
        (data["O"] / length) ** (-data["O"] / length)) + log2((data["P"] / length) ** (-data["P"] / length)) + log2(
        (data["Q"] / length) ** (-data["Q"] / length)) + log2(
        (data["R"] / length) ** (-data["R"] / length)) + log2((data["S"] / length) ** (-data["S"] / length)) + log2(
        (data["T"] / length) ** (-data["T"] / length)) + log2(
        (data["U"] / length) ** (-data["U"] / length)) + log2((data["V"] / length) ** (-data["V"] / length)) + log2(
        (data["W"] / length) ** (-data["W"] / length)) + log2((data["Y"] / length) ** (-data["Y"] / length)))
    Pit.columns = ["Pit"]

    Organism_ID = pd.DataFrame(data.Organism_ID)
    Organism_Pit = Organism_ID.join(Pit)
    Proteomes = pd.read_csv(proteomes_file, sep="\t")
    Sp_pit = Organism_Pit.merge(Proteomes, left_on="Organism_ID", right_on="Organism ID").drop_duplicates('Organism_ID',
                                                                                                          keep='last')
    Sp_pit['Kingdom'] = Sp_pit['Taxonomic lineage'].str.split(', ', expand=True)[0]
    Sp_pit['Taxonomic lineage'] = Sp_pit['Taxonomic lineage'].str.split(', ', expand=True)[1]
    Sp_pit = Sp_pit[Sp_pit.Kingdom != 'Viruses'].drop(columns=['Kingdom', 'Proteome ID', 'Organism ID']).rename(
        columns={'Taxonomic lineage': 'Phylum'}).reset_index(drop=True)

    code = []
    for row in Sp_pit['Phylum']:
        if row == "Proteobacteria":
            code.append("1")
        elif row == "Candidatus Hydrogenedentes":
            code.append("2")
        elif row == "Candidatus Abyssubacteria":
            code.append("3")
        elif row == "Spirochaetes":
            code.append("4")
        elif row == "Deferribacteres":
            code.append("5")
        elif row == "Chrysiogenetes":
            code.append("6")
        elif row == "Acidobacteria":
            code.append("7")
        elif row == "Thermodesulfobacteria":
            code.append("8")
        elif row == "Nitrospirae":
            code.append("9")
        elif row == "Nitrospinae/Tectomicrobia group":
            code.append("10")
        elif row == "Elusimicrobia":
            code.append("11")
        elif row == "Candidatus Omnitrophica":
            code.append("12")
        elif row == "Planctomycetes":
            code.append("13")
        elif row == "Chlamydiae":
            code.append("14")
        elif row == "Lentisphaerae":
            code.append("15")
        elif row == "Candidatus Aureabacteria":
            code.append("16")
        elif row == "Kiritimatiellaeota":
            code.append("17")
        elif row == "Verrucomicrobia":
            code.append("18")
        elif row == "Candidatus Latescibacteria":
            code.append("19")
        elif row == "Gemmatimonadetes":
            code.append("20")
        elif row == "Candidatus Fermentibacteria":
            code.append("21")
        elif row == "Candidatus Marinimicrobia":
            code.append("22")
        elif row == "Calditrichaeota":
            code.append("23")
        elif row == "Rhodothermaeota":
            code.append("24")
        elif row == "Balneolaeota":
            code.append("25")
        elif row == "Ignavibacteriae":
            code.append("26")
        elif row == "Candidatus Kryptonia":
            code.append("27")
        elif row == "Chlorobi":
            code.append("28")
        elif row == "Bacteroidetes":
            code.append("29")
        elif row == "Candidatus Kapabacteria":
            code.append("30")
        elif row == "Candidatus Cloacimonetes":
            code.append("31")
        elif row == "Fibrobacteres":
            code.append("32")
        elif row == "Synergistetes":
            code.append("33")
        elif row == "Fusobacteria":
            code.append("34")
        elif row == "Deinococcus-Thermus":
            code.append("35")
        elif row == "Coprothermobacterota":
            code.append("36")
        elif row == "Thermotogae":
            code.append("37")
        elif row == "Aquificae":
            code.append("38")
        elif row == "Dictyoglomi":
            code.append("39")
        elif row == "Firmicutes":
            code.append("40")
        elif row == "Tenericutes":
            code.append("41")
        elif row == "Candidatus Eremiobacteraeota":
            code.append("42")
        elif row == "Abditibacteriota":
            code.append("43")
        elif row == "Armatimonadetes":
            code.append("44")
        elif row == "Thermobaculum":
            code.append("45")
        elif row == "Chloroflexi":
            code.append("46")
        elif row == "Candidatus Dormibacteraeota":
            code.append("47")
        elif row == "Actinobacteria":
            code.append("48")
        elif row == "Cyanobacteria":
            code.append("49")
        elif row == "Candidatus Melainabacteria":
            code.append("50")
        elif row == "Candidatus Saganbacteria":
            code.append("51")
        elif row == "Candidatus Margulisbacteria":
            code.append("52")
        elif row == "candidate division KSB1" or row == "candidate division KSB3" or row == "candidate division NC10" or row == "Candidatus Aerophobetes" or row == "Candidatus Aminicenantes" or row == "Candidatus Atribacteria" or row == "Candidatus Bipolaricaulota" or row == "Candidatus Buchananbacteria" or row == "Candidatus Dadabacteria" or row == "Candidatus Delongbacteria" or row == "Candidatus Edwardsbacteria" or row == "Candidatus Eisenbacteria" or row == "Candidatus Firestonebacteria" or row == "Candidatus Goldbacteria" or row == "Candidatus Hydrothermae" or row == "Candidatus Kerfeldbacteria" or row == "Candidatus Komeilibacteria" or row == "Candidatus Lindowbacteria" or row == "Candidatus Magasanikbacteria" or row == "Candidatus Peregrinibacteria" or row == "Candidatus Poribacteria" or row == "Candidatus Riflebacteria" or row == "Candidatus Rokubacteria" or row == "Candidatus Saccharibacteria" or row == "Candidatus Schekmanbacteria" or row == "Candidatus Sumerlaeota" or row == "Candidatus Terrybacteria" or row == "Candidatus Uhrbacteria" or row == "Candidatus Wallbacteria" or row == "Candidatus Yanofskybacteria" or row == "unclassified Parcubacteria group":
            code.append("53")
        elif row == "Candidatus Altiarchaeota":
            code.append("54")
        elif row == "Candidatus Micrarchaeota":
            code.append("55")
        elif row == "Candidatus Diapherotrites":
            code.append("56")
        elif row == "Candidatus Aenigmarchaeota":
            code.append("57")
        elif row == "Candidatus Huberarchaea":
            code.append("58")
        elif row == "Nanoarchaeota":
            code.append("59")
        elif row == "Candidatus Parvarchaeota":
            code.append("60")
        elif row == "Candidatus Pacearchaeota":
            code.append("61")
        elif row == "Candidatus Woesearchaeota":
            code.append("62")
        elif row == "Euryarchaeota":
            code.append("63")
        elif row == "Candidatus Bathyarchaeota":
            code.append("64")
        elif row == "Thaumarchaeota":
            code.append("65")
        elif row == "Candidatus Korarchaeota":
            code.append("66")
        elif row == "Candidatus Verstraetearchaeota":
            code.append("67")
        elif row == "Candidatus Marsarchaeota":
            code.append("68")
        elif row == "Crenarchaeota":
            code.append("69")
        elif row == "Asgard group":
            code.append("70")
        elif row == "Amoebozoa":
            code.append("71")
        elif row == "Apusozoa":
            code.append("72")
        elif row == "Nucleariidae and Fonticula group":
            code.append("73")
        elif row == "Fungi":
            code.append("74")
        elif row == "Ichthyosporea":
            code.append("75")
        elif row == "Choanoflagellata":
            code.append("76")
        elif row == "Metazoa":
            code.append("77")
        elif row == "Heterolobosea":
            code.append("78")
        elif row == "Euglenozoa":
            code.append("79")
        elif row == "Parabasalia":
            code.append("80")
        elif row == "Diplomonadida":
            code.append("81")
        elif row == "Viridiplantae":
            code.append("82")
        elif row == "Rhodophyta":
            code.append("83")
        elif row == "Haptista":
            code.append("84")
        elif row == "Cryptophyta":
            code.append("85")
        elif row == "Rhizaria":
            code.append("86")
        elif row == "Alveolata":
            code.append("87")
        elif row == "Stramenopiles":
            code.append("88")
        else:
            code.append("0")

    Sp_pit['code'] = code
    Sp_pit.to_csv(outfile, index=False)
    endtime = time.time()
    timex = endtime - starttime
    print(
        "total time taken is {:3.0f} hours {:3.0f} minutes {:3.0f} seconds {:3.0f} milisecs".format((timex // 3600), (
                timex % 3600) // 60,
                                                                                                    (timex % 60 // 1),
                                                                                                    (timex - int(
                                                                                                        timex))))
