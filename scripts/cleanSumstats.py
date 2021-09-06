#!/usr/bin/python
# written by helena gaspar hagax8@gmail.com
import csv
import sys
import re
import gzip
#print '*******Usage:*******************************'
#print './nameofProgram.py input_gwas_sum_statistics_with_SNPID_P_INFO_A1_A2_fields output_name\n'
#print '*******Number of arguments:*****************'
#print len(sys.argv)-1, 'arguments.\n'
#print '*******Program name and argument list:******'
#print str(sys.argv)
#print ''

infile = sys.argv[1]
outfile = sys.argv[2]
countOrig = 0
countSNP = 0
addse = False
addN = False
setdaner = False
setcontrolmaf = False
fin = gzip.open(infile, 'rb')
outlog = open(outfile+".log", 'w')
# with gzip.open(infile, 'rb') as fin:
try:
    headline = next(fin).decode('ascii').replace(",", " ").strip().split()
    headers = [x.upper() for x in headline]
    ncols = len(headers)
    outlog.write("\n")
    outlog.write("********GWAS Header:*********\n")
    for item in headers:
        outlog.write(" %s " % item)
    outlog.write("\n")
    isInfo = False
    Infodic = ["INFO", "IMPQUALITY", "INFO.PLINK", "INFO_UKBB"]
    # Pdic=["PVALUE","P_VALUE","P","PVAL","P.VAL","GWAS_P","P-VALUE","P-VAL","FREQUENTIST_ADD_PVALUE","P.VALUE","P_VAL","SCAN-P","P.LMM","META.PVAL","P_RAN","P.ADD","P_ADJ","P_DGC","GC_PVALUE","GC.PVALUE","PVAL.GC","P-VALUE_GC","P_GC","P.2GC","P_BOLT_LMM"]
    Pdic = ["PVALUE", "P_VALUE", "P", "PVAL", "P.VAL", "GWAS_P",
            "P-VALUE", "P-VAL", "FREQUENTIST_ADD_PVALUE",
            "P.VALUE", "P_VAL", "SCAN-P", "P.LMM", "META.PVAL",
            "P_RAN", "P.ADD", "P_BOLT_LMM"]
    SNPdic = ["SNP", "SNPID", "MARKERNAME", "MARKER_NAME",
              "SNPTESTID", "ID_DBSNP49", "RSID", "ID", "RS_NUMBER",
              "MARKER", "RS", "RSNUMBER", "RS_NUMBERS", "SNP.NAME",
              "SNP ID", "SNP_ID", "LOCATIONALID", "ASSAY_NAME"]
    A1dic = ["A1", "A1_EFFECT", "REF", "EFFECT_ALLELE", "RISK_ALLELE",
             "EFFECTALLELE", "EFFECT_ALL", "REFERENCE_ALLELE",
             "REF_ALLELE", "REFERENCEALLELE", "EA", "ALLELE_1",
             "INC_ALLELE", "ALLELE1", "A", "A_1", "CODED_ALLELE", "TESTED_ALLELE"]
    A2dic = ["A2", "ALT", "A2_OTHER", "NONREF_ALLELE", "NEFFECT_ALLELE",
             "NEFFECTALLELE", "NONEFFECT_ALLELE", "OA", "OTHER_ALL",
             "ALLELE_2", "OTHER_ALLELE", "OTHERALLELE", "NEA",
             "NON_EFFECT_ALLELE", "NONEFFECTALLELE", "DEC_ALLELE",
             "NEA", "ALLELE2", "ALLELE0", "ALLELE_0", "ALT_ALLELE",
             "A_0", "NONCODED_ALLELE"]
# this is a maf dic including ref pop if someone wants to use it instead of MAFdic
#	MAFdic=["FREQ_HAPMAP","HAPMAP_FREQ","CEUAF","CEUMAF","EUR_MAF","MAF_HAPMAP_REL27","FREQ.ALLELE1.HAPMAPCEU","ALLELE1_FREQ_HAPMAPCEU","FREQ_EUROPEAN_1000GENOMES","FRQ_ALT.1KGP_ASN_N286","ALLELE1_FREQ_HAPMAPCEU"]#notclear: FREQ_A1, A1FREQ, FREQ1
    MAFdic = ["FRQ", "MAF", "EAF", "EFFECT_ALLELE_FREQ", "FREQ",
              "FREQ_U", "A1FREQ", "MA_FREQ", "MAF_NW", "FREQ_A1",
              "A1FREQ", "FREQ1", "AF", "CODED_ALLELE_FREQUENCY",
              "FREQ_TESTED_ALLELE_IN_HRS", "EAF_HRC"]
    FRQAdic = ["FRQ_A"]
    FRQUdic = ["FRQ_U"]
    Betadic = ["B", "BETA", "EFFECT_BETA", "EFFECT", "GWAS_BETA",
               "EFFECT_A1", "EFFECTA1", "EFFECT_NW"]
    Zdic = ["Z-SCORE", "ZSCORE", "BETAZSCALE", "Z"]
    ORdic = ["OR", "ODDS-RATIO", "ODDS_RATIO", "ODDSRATIO",
             "OR(MINALLELE)", "OR.LOGISTIC", "OR_RAN", "OR(A1)"]
    SEdic = ["SE", "STDERR", "STD", "STDER", "STANDARD_ERROR", "OR_SE"
             "STANDARDERROR", "STDERR_NW", "META.SE", "SE_DGC", "SE.2GC"]
    Ndic = ["TOTALN", "N_TOTAL", "N_SAMPLES", "N_ANALYZED",
            "NSAMPLES", "SAMPLESIZE", "SAMPLE_SIZE", "TOTAL_SAMPLE_SIZE",
            "TOTALSAMPLESIZE", "N", "EUROPEAN_N"]
    Chrdic = ["CHR", "CHROMOSOME", "CHR_BUILD37", "CHR_BUILD36", "CHR_B37",
              "CHR_B38", "CHR_ID", "SCAFFOLD", "HG18CHR", "CHR.HG18",
              "CHR_HG18", "CHR_BP_HG19B37", "HG19CHRC"]
    BPdic = ["BP", "POSITION", "LOCATION", "POS", "POS_B37", "POS_BUILD37",
             "CHR_POSITION", "BP_HG19B37", "POS_B36", "POS_BUILD36",
             "PHYSPOS", "POS.HG18", "POS_HG18", "BP_HG19", "BP.GRCH37",
             "GENPOS", "POSITION(HG19)", "POS(B37)"]
    NCASdic = ["NCAS", "N_CAS", "NCASES", "N_CASES",
               "N_CASE", "NCASE", "CASES", "NCA"]
    NCONdic = ["NCON", "N_CON", "NCONTROLS", "N_CONTROLS",
               "NCONTROL", "N_CONTROL", "CONTROLS", "NCO"]
    namedic = []
    indexdic = []
    defP = list(set(headers).intersection(Pdic))
    defSNP = list(set(headers).intersection(SNPdic))
    defA1 = list(set(headers).intersection(A1dic))
    defA2 = list(set(headers).intersection(A2dic))
    defMAF = list(set(headers).intersection(MAFdic))
    defInfo = list(set(headers).intersection(Infodic))
    defOR = list(set(headers).intersection(ORdic))
    defBeta = list(set(headers).intersection(Betadic))
    defSE = list(set(headers).intersection(SEdic))
    defN = list(set(headers).intersection(Ndic))
    defChr = list(set(headers).intersection(Chrdic))
    defBP = list(set(headers).intersection(BPdic))
    defNcas = list(set(headers).intersection(NCASdic))
    defNcon = list(set(headers).intersection(NCONdic))
    defZ = list(set(headers).intersection(Zdic))
    defFRQA = list(set(headers).intersection(FRQAdic))
    defFRQU = list(set(headers).intersection(FRQUdic))
    outlog.write("\n")
# check fields
    if not defChr:
        outlog.write("No Chromosome field\n")
    else:
        outlog.write("Chromosome <=> %s\n" % defChr[0])
        chrindex = headers.index(defChr[0])
        namedic.append("CHR")
        indexdic.append(chrindex)

    if not defSNP:
        outlog.write("No SNPID field\n")
        quit()
    else:
        outlog.write("SNP <=> %s\n" % defSNP[0])
        snpindex = headers.index(defSNP[0])

    if not defBP:
        outlog.write("no BP field\n")
    else:
        outlog.write("BP <=> %s\n" % defBP[0])
        bpindex = headers.index(defBP[0])
        namedic.append("ORIGBP")
        indexdic.append(bpindex)

    if not defA1:
        outlog.write("No A1 field\n")
    else:
        outlog.write("A1 <=> %s\n" % defA1[0])
        a1index = headers.index(defA1[0])
        namedic.append("A1")
        indexdic.append(a1index)

    if not defA2:
        outlog.write("No A2 field\n")
    else:
        outlog.write("A2 <=> %s\n" % defA2[0])
        a2index = headers.index(defA2[0])
        namedic.append("A2")
        indexdic.append(a2index)

    if not defMAF and not defFRQA and not defFRQU:
        try:
            frqa = filter(lambda x: x.startswith('FRQ_A_'), headers)[0]
            frqu = filter(lambda x: x.startswith('FRQ_U_'), headers)[0]
            if frqu and frqa:
                outlog.write(
                    "Found daner-like MAF fields: %s and %s\n" % (frqa, frqu))
                Ncas = float(frqa[6:])
                Ncon = float(frqu[6:])
                N = Ncas+Ncon
                outlog.write("Ncas, Ncon, Ntot = %s, %s, %s\n" %
                             (Ncas, Ncon, N))
                setdaner = True
                setcontrolmaf = True
                mafindex = headers.index(frqu)
                namedic.append("FREQ")
                indexdic.append(mafindex)
                namedic.append("FREQ_CASES")
                mafcaseindex = headers.index(frqa)
                indexdic.append(mafcaseindex)
        except:
            outlog.write("No daner-like field found")
    elif defFRQA and defFRQU:
        outlog.write("Found case/control MAF fields: %s and %s\n" %
                     (defFRQA[0], defFRQU[0]))
        namedic.append("FREQ")
        mafindex = headers.index(defFRQU[0])
        indexdic.append(mafindex)
        namedic.append("FREQ_CASES")
        mafcaseindex = headers.index(defFRQA[0])
        indexdic.append(mafcaseindex)
        setcontrolmaf = True
    else:
        outlog.write("Found MAF field: %s\n" % defMAF[0])
        mafindex = headers.index(defMAF[0])
        namedic.append("FREQ")
        indexdic.append(mafindex)

    if not defP:
        outlog.write("No P field\n")
    else:
        outlog.write("P <=> %s\n" % defP[0])
        pindex = headers.index(defP[0])
        namedic.append("P")
        indexdic.append(pindex)

    if not defBeta:
        outlog.write("No Beta field\n")
        if not defZ:
            outlog.write("No Zscore field either\n")
        else:
            outlog.write(
                "Beta (here zscore, SE will be set to 1) <=> %s\n" % defZ[0])
            zindex = headers.index(defZ[0])
            namedic.append("BETA")
            indexdic.append(zindex)
            addse = True
    else:
        outlog.write("Beta <=> %s\n" % defBeta[0])
        betaindex = headers.index(defBeta[0])
        namedic.append("BETA")
        indexdic.append(betaindex)


    if not defOR:
        outlog.write("No OR field\n")
    elif not defBeta and not defZ:
        outlog.write("OR <=> %s\n" % defOR[0])
        orindex = headers.index(defOR[0])
        namedic.append("OR")
        indexdic.append(orindex)

    if addse == True:
        outlog.write("SE set to 1.\n")
    elif not defSE:
        outlog.write("No SE field.\n")
    else:
        outlog.write("SE <=> %s\n" % defSE[0])
        seindex = headers.index(defSE[0])
        namedic.append("SE")
        indexdic.append(seindex)

    if not defNcon:
        outlog.write("No Ncon field\n")
    else:
        setdaner = False
        outlog.write("Ncon <=> %s\n" % defNcon[0])
        nconindex = headers.index(defNcon[0])
        namedic.append("Ncon")
        indexdic.append(nconindex)

    if not defNcas:
        outlog.write("No Ncas field\n")
    else:
        setdaner = False
        outlog.write("Ncas <=> %s\n" % defNcas[0])
        ncasindex = headers.index(defNcas[0])
        namedic.append("Ncas")
        indexdic.append(ncasindex)

    if not defN:
        outlog.write("No N field\n")
        if defNcon and defNcas:
            outlog.write("Will be inferring N field from Ncon and Ncas\n")
            addN = True
    else:
        outlog.write("N <=> %s\n" % defN[0])
        nindex = headers.index(defN[0])
        namedic.append("N")
        indexdic.append(nindex)

    if not defInfo:
        outlog.write("No INFO field\n")
    else:
        outlog.write("INFO <=> %s\n" % defInfo[0])
        infoindex = headers.index(defInfo[0])
        namedic.append("INFO")
        indexdic.append(infoindex)

    outlog.write("\n")
    if defSNP and defP and (defBeta or defZ or defOR) and defA1 and defA2:
        outlog.write("ENOUGH_FOR_PRS\n")
    if defSNP and defP:
        outlog.write("ENOUGH_FOR_PATHWAY_ANALYSIS\n")
    outlog.write("\n")

    outlog.write("Variant filtering/cleaning. It could take 2-3 minutes.\n")
    outlog.close()
    mylogstring = ""
    out = gzip.open(outfile+".gz", 'wb')
    try:
        for line in fin:
            row = line.decode('ascii').replace(",", " ").strip().split()
            countOrig += 1
            if len(row) != ncols:
                mylogstring += "ERROR: incorrect number of fields at " + \
                    str(countOrig) + "\n"
            else:
                if (countOrig == 1):
                    out.write(b'SNP')
                    for i in namedic:
                        out.write(("\t%s" % (i)).encode())
                    if addse:
                        out.write(b"\tSE")
                    if addN:
                        out.write(b"\tN")
                    if setdaner:
                        out.write(b"\t%s\t%s\t%s" % ("Ncas", "Ncon", "N"))
                    out.write(b"\n")
                try:
                    cleanRSID = re.sub('^X:', '23:', re.sub(
                        '^chr', '', re.sub('_', ':', row[snpindex])))
                except:
                    outlog.write("no id at row %s\n" % countOrig)
                    cleanRSID = b"NORSID"
                    pass
                countSNP += 1
                out.write(cleanRSID.encode())
                for i in indexdic:
                    out.write(("\t%s" % row[i].upper()).encode())
                if addse:
                    out.write(b"\t1")
                if addN:
                    out.write(("\t%f" %
                              (float(row[ncasindex])+float(row[nconindex]))).encode())
                if setdaner:
                    out.write(("\t%s\t%s\t%s" % (Ncas, Ncon, N)).encode())
                out.write(b"\n")
    finally:
        out.close()
        outlog = open(outfile+".log", 'a')
        outlog.write(mylogstring)
        outlog.write("Initial number of variants: %s\n" % countOrig)
        outlog.write("Final number of variants: %s\n" % countSNP)
        outlog.write("Removed: %s\n" % (int(countOrig)-int(countSNP)))
finally:
    fin.close()
    outlog.close()
