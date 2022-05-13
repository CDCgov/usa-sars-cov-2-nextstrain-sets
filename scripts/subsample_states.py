#!/usr/bin/env python

### Written by Jason Caravas, mqi4 
### Create subsamples for ncbi state data

# from impala.dbapi import connect
import os, sys, re, argparse, csv, lzma, pyodbc
import datetime as dt
import pandas as pd
import numpy as np

##  Download data from hadoop database
def loadUSAData(usaData, outDir):
    ##  Select fields from USA isolates with associated state data and updated pangolin calls

    df = pd.read_csv(usaData, sep = "\t")
    df.columns = ["accession", "organism", "strain", "isolate", "isolation_source", "host", "country", "collection_date", "pub_date", "update_date", "bioprojects", "biosample", "completeness", "lineage", "clade", "seq_length", "nt_id"]

    ##  Remove duplicate entries from query
    df = df.drop_duplicates()

    outFile = os.path.join (outDir, "1_raw_metadata.tsv")
    df.to_csv(outFile, sep="\t")
    return df

##  Download reference strains from hadoop database
def loadRefData(refData, outDir):
    ##  Select fields from Wuhan reference strains with associated state data and updated pangolin calls

    df = pd.read_csv(refData, sep = "\t")
    df.columns = ["accession", "organism", "strain", "isolate", "isolation_source", "host", "country", "collection_date", "pub_date", "update_date", "bioprojects", "biosample", "completeness", "lineage", "clade", "seq_length", "nt_id"]

    ##  Note: Wuhan correction code also present in downloadSequences function
    ##  Add common name to WHO1 strain name (in "isolate" field for now).  It is null in ncbi ingest
    df.loc[df.accession == "LR757998", "isolate"] = "Wuhan/WH01/2019"
    ##  Rename Wuhan-Hu-1 to match name hardcoded in Nextstrain snakemake files
    df.loc[df.accession == "NC_045512", "isolate"] = "Wuhan/Hu-1/2019"
    
    ##  Set Wuhan-Hu-1 to default date of 2019-12-01 (only month specified)
    df.loc[df.accession == "NC_045512", "collection_date"] = "2019-12-01"
    
    
    ##  Convert dataframe to nextstrain fields
    df["normalized_country"] = "Hubei"
    df.rename( mapper =
                        {   "accession"             :   "genbank_accession",
                            "organism"              :   "virus",
                            "country"               :   "original_country",
                            "normalized_country"    :   "division",
                            "collection_date"       :   "date",
                            "pub_date"              :   "date_submitted",
                            "lineage"               :   "pango_lineage",
                            "seq_length"            :   "length",
                            "clade"                 :   "Nextstrain_clade"
                        },
                        inplace = True,
                        axis = 1
                      )
    ## Fill in nextstrain fields to best of my ability
    ## Static fields:
    df["virus"] = "ncov"
    df["gisaid_epi_isl"] = "?"
    df["region"] = "Asia"
    df["location"] = np.nan
    df["country"] = "China"
    df["region_exposure"] = "Asia"
    df["country_exposure"] = "China"
    df["segment"] = "genome"
    df["age"] = "?"
    df["sex"] = "?"
    df["GISAID_clade"] = np.nan
    df["originating_lab"] = "?"
    df["submitting_lab"] = "?"
    df["authors"] = "?"
    df["url"] = "?"
    df["title"] = "?"
    df["paper_url"] ="?"
    df["purpose_of_sequencing"] ="?"
    
    ##  Dependent Fields:
    df["division_exposure"] = np.nan
    df["division_exposure"] = df["division_exposure"].fillna(df["division"])
    
    ##  Overwrite "strain" value with "isolate" value
    df["strain"] = df["isolate"]

    outFile = os.path.join (outDir, "1_Wuhan_metadata.tsv")
    df.to_csv(outFile, sep="\t")
    return df
    
def downloadSequences (nextstrainDf, outDir):
    idList = nextstrainDf["genbank_accession"].tolist()
    
    conn = pyodbc.connect(connInfo, autocommit=True)
    ##  Select fields from nextstrainDf with associated state data and updated pangolin calls
    sequenceDf = pd.DataFrame(columns = ["accession", "strain", "isolate", "sequence"])
    done = False
    maxQuery = 9999
    end = maxQuery
    i = 0
    while not done:
        if i + maxQuery < len (idList):
            queryList = idList[i:i + maxQuery]
            queryListString = ','.join(map("'{0}'".format, queryList))
            selectStatement =   "SELECT ncbi_ingest.accession, ncbi_ingest.strain, ncbi_ingest.isolate, ncbi_ingest.sequence FROM ncbi_ingest WHERE ncbi_ingest.accession IN (" + queryListString + ")"
            sequenceDf = sequenceDf.append ( pd.read_sql(selectStatement, conn) )
            i += maxQuery
        else:
            queryList = idList[i:]
            queryListString = ','.join(map("'{0}'".format, queryList))
            selectStatement =   "SELECT ncbi_ingest.accession, ncbi_ingest.strain, ncbi_ingest.isolate, ncbi_ingest.sequence FROM ncbi_ingest WHERE ncbi_ingest.accession IN (" + queryListString + ")"
            sequenceDf = sequenceDf.append ( pd.read_sql(selectStatement, conn) )
            done = True

    ##  Note: Wuhan correction code also present in downloadRefData function
    ##  Add common name to WHO1 strain name (in "isolate" field for now).  It is null in ncbi ingest
    sequenceDf.loc[sequenceDf.accession == "LR757998", "isolate"] = "Wuhan/WH01/2019"
    ##  Rename Wuhan-Hu-1 to match name hardcoded in Nextstrain snakemake files
    sequenceDf.loc[sequenceDf.accession == "NC_045512", "isolate"] = "Wuhan/Hu-1/2019"
    
    ##  Overwrite "strain" value with "isolate" value
    sequenceDf["strain"] = sequenceDf["isolate"]
        
    outFile = os.path.join (outDir, "sequences.tsv")
    sequenceDf.to_csv(outFile, sep="\t")
    conn.close()
    return sequenceDf

##  Eliminate all entries missing essential metadata and below minimum sequence length
def filterDataframeMissing(rawDf, minLength, outDir):
    ##  Notes on non-required fields
    ##  "strain": isolate is better than strain and strain value will be overwritten by isolate
    ##          in prepareFiles subroutine
    ##  "isolation_source:  Not necessary for nextstrain.
    ##  "update_date:   Not necessary for nextstrain
    ##  "bioprojects":  Not necessary for nextstrain
    ##  "biosample":    Not necessary for nextstrain
    ##  "completeness": Not necessary for nextstrain.  Not very informative
    
    ##  Notes on "lineage":
    ##  The NCBI pangolin lineage is not used and instead the newest call from the 
    ##  pangolin table is used (see select statement)
    
    requiredFields = ["accession", "organism", "isolate", "host", "country", "collection_date", "pub_date", "lineage", "seq_length", "nt_id"]
    
    metadataCompleteDf = rawDf.copy()
    
    ##  Simple method to do length filter here is to unassign seq_length values
    ##  below minLength threshold and then let the Df.dropna() function eliminate 
    ##  them in the next step
    ##  Unassign seq_length < minLength
    ##  https://www.activestate.com/resources/quick-reads/how-to-apply-functions-in-pandas/
    metadataCompleteDf["seq_length"] = metadataCompleteDf["seq_length"].apply(lambda val: np.nan if val < minLength else val)
    
    ##  Flag entries with misformatted date as null for removal
    metadataCompleteDf["collection_date"].replace([r'^(?!.*\d{4}-\d{2}-\d{2}).*$'], np.nan, regex = True, inplace = True)
    
    metadataCompleteDf = metadataCompleteDf.dropna(subset = requiredFields, how = 'any')
    
    outFile = os.path.join (outDir, "2_complete_metadata.tsv")
    metadataCompleteDf.to_csv(outFile, sep="\t")
    
    return metadataCompleteDf
    
##  Normalize field name for state data
def normalizeStateData (metadataCompleteDf, normFile, outDir):
    stateAbbreviationLookup = []
    
    normalizedDf = metadataCompleteDf.copy()
    ##  Read in state/terr abbreviations
    with open(normFile, newline='') as csvfile:
        csvReader = csv.reader(csvfile,delimiter='\t')
        for row in csvReader:
            stateAbbreviationLookup.append(row)

    ##  Precompile splitter regex
    splitRegEx = re.compile("\s*[,\/]\s*")
    
    ##  Apply lambda to normalize state (country) field
    ##  https://www.activestate.com/resources/quick-reads/how-to-apply-functions-in-pandas/
    normalizedDf["normalized_country"] = normalizedDf.apply(lambda row: normalizeStateFieldLambda(row["country"], stateAbbreviationLookup, splitRegEx), axis = 1)

    ##  Eliminate entries with missing state data
    normalizedDf = normalizedDf.dropna(subset = ["normalized_country"], how = 'any')
    
    outFile = os.path.join (outDir, "3_normalized_metadata.tsv")
    normalizedDf.to_csv(outFile, sep="\t")
    return normalizedDf

##  Lambda function for pandas dataframe to add a normalized ormalize state (country) field
def normalizeStateFieldLambda (rawCountry, stateAbbreviationLookup, splitRegEx):
    ##  Remove "USA: " prefix
    matches = re.search("USA:\s*(.+)", rawCountry)
    state = matches.group(1)
    
    
    ##  Skip unknowns
    if state == "Unknown":
        print (f"{rawCountry} could not be matched to state/territory.  Excluding.", file=sys.stderr)
        return np.nan
    
    ##  Split entries with multiple levels of specificity
    matches = splitRegEx.split(state)
                
    normalizedState = ""
    for field in matches:
        for entry in stateAbbreviationLookup:
            if field.lstrip() in entry:
                normalizedState = entry[0]
                return normalizedState
    if not normalizedState:
        print (f"{rawCountry} could not be matched to state/territory.  Excluding.", file=sys.stderr)
        return np.nan
        
##  Eliminate identical sequences within same state
def removeNonUniqueSequences(normalizedDf, outDir):
    uniqueDf = normalizedDf.copy()
    
    uniqueDf = uniqueDf.drop_duplicates(subset = ["normalized_country", "nt_id"], keep = "first") 
    
    outFile = os.path.join (outDir, "4_unique_metadata.tsv")
    uniqueDf.to_csv(outFile, sep="\t")
    return uniqueDf
    
##  Create full state dumps
def fullSample(state, normalizedDf, refDf, outDir):
    subsampleLists = ""
   
    ##  Create state specific dataframe
    stateDf = normalizedDf.loc[normalizedDf["normalized_country"] == state].copy()
    
    ##  Prepare output
    sampleCount = len(stateDf.index)
    if sampleCount == 0:
        print(f"{sampleCount} samples found for {state}.  Terminating with error.")
        sys.exit(1)
    
    stateFileRoot=state.replace(" ", "")
    stateDir = os.path.join(outDir, stateFileRoot)
    if not os.path.exists(stateDir):
        os.makedirs(stateDir)
    sampleName = stateFileRoot + "-full"
        
    prepareFiles (stateDf, refDf, sampleName, state, outDir)
        
##  Create random subsample dataframes for each state
def subsampleUnweighted(state, uniqueDf, refDf, sampleSizes, outDir):
    subsampleLists = ""

    stateFileRoot=state.replace(" ", "")
    stateDir = os.path.join(outDir, stateFileRoot)
    if not os.path.exists(stateDir):
        os.makedirs(stateDir)
    
    ##  Create state specific dataframe
    stateDf = uniqueDf.loc[uniqueDf["normalized_country"] == state].copy()
    
    ##  Shuffle dataframe in-place using .sample()
    stateDf = stateDf.sample(frac=1).reset_index(drop=True)

    sampleSizes.sort()
   
    # Generate nested subsamples from randomized list
    for size in sampleSizes:
        stateDfCopy = stateDf.copy()
        sampleCount = len(stateDfCopy.index)
        sampleName = stateFileRoot + "-" + str(size) + "-unweighted"
        if sampleCount > size:
            subsampleDf = stateDfCopy[:size].copy()
            prepareFiles (subsampleDf, refDf, sampleName, state, outDir)
        else:
            prepareFiles (stateDfCopy, refDf, sampleName, state, outDir)
            
##  Create random subsample dataframes for each state weighted towards recent samples
##      Scaling exponent controls the strength of weighting effect.
##      Higher scaling factor = stronger skew towards recent samples.
##      Isolates will be binned based on the the number of <period> day periods prior to
##      the current date the sample was collected, starting with bin "1" (not "0").
##      Relative weighting of entire date bin is
##          (1/bin) exp(scalingExponent)
##      with individual isolates weighted as
##          bin weight / count of isolates in bin
##      Absolute weighting is determined by summing isolate weights and 
##      dividing each weight by that sum, normalizing total weight to 1.
##      Note that entire date bins are weighted together, so weighted sample should be
##      independent of number of sequences generated within the binned time window.
def subsampleWeighted(state, uniqueDf, refDf, sampleSizes, scalingExponent, period, current_date, outDir):
    subsampleLists = ""
    stateFileRoot=state.replace(" ", "")
    stateDir = os.path.join(outDir, stateFileRoot)
    if not os.path.exists(stateDir):
        os.makedirs(stateDir)
    
    ##  Create state specific dataframe
    stateDf = uniqueDf.loc[uniqueDf["normalized_country"] == state].copy()
    
    ##  Assign bins
    stateDf["bin"] = stateDf.apply(lambda row: assignBinsLambda(row["collection_date"], current_date, period), axis = 1)
    
    ##  Assign absolute weight based on scaling factors and number of items in bin
    valueTotal = 0
    for bin in stateDf["bin"].unique():
        relativeBinWeight = (1 / bin) ** scalingExponent
        valueTotal += relativeBinWeight
    
    stateDf["absolute_weight"] = np.nan
    
    for bin in stateDf["bin"].unique():
        count = stateDf[stateDf["bin"] == bin]["accession"].count()
        relativeBinWeight = (1 / bin) ** scalingExponent
        absoluteIsolateWeight = (relativeBinWeight / count) / valueTotal
        stateDf.loc[stateDf.bin == bin, "absolute_weight"] = absoluteIsolateWeight
       
    sampleSizes.sort()
   
    ##  Generate weighted subsamples from list
    ##  Weighted subsamples are not nested (smaller samples are not a subset of larger samples)
    for size in sampleSizes:
        stateDfCopy = stateDf.copy()
        subsampleCount = size
        sampleName = stateFileRoot + "-" + str(subsampleCount) + "-weighted"
        if (len(stateDfCopy.index) > size ):            
            subsampleDf = stateDfCopy.sample(n = subsampleCount, weights = stateDfCopy.absolute_weight, axis = 0)
            prepareFiles (subsampleDf, refDf, sampleName, state, outDir)
        else:
            prepareFiles (stateDfCopy, refDf, sampleName, state, outDir)
       
##  Lambda to assign bin.  See "subsampleWeighted" subroutine
##  for description
def assignBinsLambda (isoDateString, currentDate, period):
    isoDate = dt.datetime.strptime(isoDateString, "%Y-%m-%d")
    timedelta = currentDate-isoDate
    days = timedelta.days
    bin = (int ( round ( days / period ) ) ) + 1
    return bin

##  Create subsample data files
def prepareFiles (subsampleDf, refDf, sampleName, state, outDir):
    stateFileRoot=state.replace(" ", "")
    metadataFile = os.path.join (outDir, stateFileRoot, sampleName + "-metadata.tsv.xz")
    fastaFile = os.path.join (outDir, stateFileRoot, sampleName + "-sequences.fasta.xz")
    
    ##  Convert dataframe to nextstrain fields
    subsampleDf.rename( mapper =
                        {   "accession"             :   "genbank_accession",
                            "organism"              :   "virus",
                            "country"               :   "original_country",
                            "normalized_country"    :   "division",
                            "collection_date"       :   "date",
                            "pub_date"              :   "date_submitted",
                            "lineage"               :   "pango_lineage",
                            "seq_length"            :   "length",
                            "clade"                 :   "Nextstrain_clade"
                        },
                        inplace = True,
                        axis = 1
                      )
    ## Fill in nextstrain fields to best of my ability
    ## Static fields:
    subsampleDf["virus"] = "ncov"
    subsampleDf["gisaid_epi_isl"] = "?"
    subsampleDf["region"] = "North America"
    subsampleDf["location"] = np.nan
    subsampleDf["country"] = "USA"
    subsampleDf["region_exposure"] = "North America"
    subsampleDf["country_exposure"] = "USA"
    subsampleDf["segment"] = "genome"
    subsampleDf["age"] = "?"
    subsampleDf["sex"] = "?"
    subsampleDf["GISAID_clade"] = np.nan
    subsampleDf["originating_lab"] = "?"
    subsampleDf["submitting_lab"] = "?"
    subsampleDf["authors"] = "?"
    subsampleDf["url"] = "?"
    subsampleDf["title"] = "?"
    subsampleDf["paper_url"] ="?"
    subsampleDf["purpose_of_sequencing"] ="?"
    
    ##  Dependent Fields:
    ##  There may be a more efficient way to do this
    subsampleDf["division_exposure"] = np.nan
    subsampleDf["division_exposure"] = subsampleDf["division_exposure"].fillna(subsampleDf["division"])
    
    ##  Overwrite "strain" value with "isolate" value
    ##  Nextstrain is dropping "SARS-CoV-2/" prefix when it combines metadata
    ##  This is causing it to fail matching fasta sequence data to metadata
    ##  Remove "SARS-CoV-2/" from fasta file and remove from metadata so files agree for user
    subsampleDf["strain"] = subsampleDf["isolate"].apply (lambda val: re.sub("^SARS-CoV-2\/", "", val))
    
    ##  Add reference metadata
    combinedDf = subsampleDf.copy()
    
    ##  Create new dataframe with nextstrain metadata columns in print order
    nextstrainDf = combinedDf[["strain", "virus", "gisaid_epi_isl", "genbank_accession", "date", "region", "country", "division", "location", "segment", "length", "host", "age", "sex", "Nextstrain_clade", "pango_lineage", "GISAID_clade", "originating_lab", "submitting_lab", "authors", "url", "title", "paper_url", "date_submitted", "purpose_of_sequencing"]]
    
    nextstrainDf.to_csv(metadataFile, sep = "\t", index = False, compression = "xz")
    
    sequenceDf = downloadSequences(nextstrainDf, outDir)
    
    ##  Apply lambda to write Fasta
    ##  https://www.activestate.com/resources/quick-reads/how-to-apply-functions-in-pandas/
    with lzma.open(fastaFile, "w") as fastaFileFh:
        sequenceDf.apply(lambda row: writeFastaLambda(row["strain"], row["sequence"], fastaFileFh), axis = 1)
    

def writeFastaLambda (strain, sequence, fh):
    ##  Nextstrain is dropping "SARS-CoV-2/" prefix when it combines metadata
    ##  This is causing it to fail matching fasta sequence data to metadata
    ##  Remove "SARS-CoV-2/" from fasta file
    strain = re.sub("^SARS-CoV-2\/", "", strain)
    isoString = ">" + strain + "\n" + sequence + "\n"
    bString = isoString.encode()
    fh.write(bString)

def main(arguments):

    parser = argparse.ArgumentParser(description = "Extract and normalize NCBI state/territory data and metadata")
    parser.add_argument("-s", "--state", type = str, help = "State or territory name", required=True)
    parser.add_argument("-u", "--usaData", type = str, help = "USA metadata file", required=True)
    parser.add_argument("-r", "--refData", type = str, help = "Wuhan reference data file", required=False)
    parser.add_argument("-n", "--normFile", type = str, nargs = '?', help = "File to be used to normalize country/state names", default = "state-terr_lookup.tsv")
    parser.add_argument("-m", "--minLength", type = int, nargs = '?', help = "Minimum genome length", default = 29000)
    parser.add_argument("-o", "--output", type = str, help="Output dir", required=True)

    
    args = parser.parse_args()
    

    ##  Subsamples cannot be greater than 10000 due to query construction limits on db
    sampleSizes = [1000, 2000, 4000]
    scalingExponent = 0.8
    period = 7
    
    ##  User must specify connection information
    connInfo = ""
    
    ##  Get current date once at start of execution so that date changes during
    ##  run time do not impact date based analysis
    currentDate = dt.datetime.today()
    
    ## Create output dir and change into it
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        
    ##  Download reference data from hadoop database
    ##  NOTE: This metadata is converted to the same format as produced by "produceFiles"
    ##      subroutine for ease of use downstream.  It is incompatible with the other
    ##      dataframe processing subroutines
    
    ## Make refDf optional.  Generate empty df placeholder and fill if refData provided
    refDf = pd.DataFrame(columns=["accession", "organism", "strain", "isolate", "isolation_source", "host", "country", "collection_date", "pub_date", "update_date", "bioprojects", "biosample", "completeness", "lineage", "clade", "seq_length", "nt_id"])

    if args.refData:
        refDf = loadRefData (args.refData, args.output)
    
    ##  Download USA data from hadoop database
    rawDf = loadUSAData(args.usaData, args.output)

    ##  Eliminate all entries missing essential metadata and below minimum sequence length
    metadataCompleteDf = filterDataframeMissing(rawDf, args.minLength, args.output)
    ##  Normalize field name for state data
    normalizedDf = normalizeStateData(metadataCompleteDf, args.normFile, args.output)
    
    
    
    ##  Complete state data dumps
    fullSample (args.state, normalizedDf, refDf, args.output)
    
    ##  Eliminate identical sequences within same state
    uniqueDf = removeNonUniqueSequences(normalizedDf, args.output)
    
    ## Create random subsamples
    subsampleUnweighted (args.state, uniqueDf, refDf, sampleSizes, args.output)
    subsampleWeighted (args.state, uniqueDf, refDf, sampleSizes, scalingExponent, period, currentDate, args.output)
    
if __name__=='__main__':
        main(sys.argv[1:])
