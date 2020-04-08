##########################################################################################################################################
##########################################################################################################################################
#                    																													 #
# Author: Luca Santini                    	    																						 #
# e-mail: luca.santini.eco@gmail.com		            																				 #
# Title: Game of Tenure: The role of "hidden" citations on researchers' ranking in ecology                                               #
# Journal: XXXXXXX              																							 #
#                                      																									 #
##########################################################################################################################################
##########################################################################################################################################

#Conceptual description of the algorithm
#The simulation is initialized with the 40 journals and the number of authors in the sample (nAu). 
#The simulation lasts 10 years, all authors start with no citations and the journals with no impact factor

#Every year the simulation:
#1.	Assigns to each author the number of papers published sampling from a Poisson distribution centered to the empirical PubRate estimated for the author.
#2.	Generates a corresponding number of papers per author and assigns them to a category. 
#The category is sampled proportionally to the PropCategories of the author. 
#3.	Assigns each paper to a journal, sampling from the frequency of journals published in each category (nPJ). 
#In order to assign main text citations, we also assign an ID to each paper that corresponds to a paper sampled from WOS for the same journal.
#4.	Assigns to each paper nRefSM sampling from the empirical distributions of the same journal.
#5.	Calculates the total number of supplementary citations to distribute to papers from the total number of supplementary references.
#6. For each journal and for present and past years (from year 1 to present):
#	a) Assigns the new main text citations to papers every year corresponding to the number of citations received from sampled papers published in the same journal. 
#	   Citations are used as means of a Poisson distribution for the citations received every year, to reduce the determinism of citations received by each paper per year in our sample.
#	b) Calculates the total supplementary citations to distribute to the papers of the journal, we multiply the total number of references in SM of all papers published in the year to pY of the respective year.
#	c) To correct for the different numbers of papers published per journal, for each journal the overall citations are divided by the proportion of papers published in the simulation over the total published in 2017. 
#	   This enables us to estimate the citations coming from other papers in the same journal that are not simulated.
#	d) For each journal, distributes the supplementary citations to the other journals using the proportion of supplementary references of the journal to the others (citSupMat).
#	e) For each receiving journal, the supplementary citations received are multiplied by the proportion of the papers in the simulation over the total published in 2017. 
#	   The first year supplementary citations are distributed randomly across papers in the journal, but from the second year the citations received by each paper are sampled from a Poisson distribution centered to the citations received the previous year + 1, thus ensuring that papers that did not get any citation the previous year still have a chance to be cited. The sum of the citations for the journal is then rescaled and rounded to the total citations the journal receives.
#	f) At the end of the year, calculates the number of total citations per author, and the H-Index per author. 
#      From the third year of the simulation, every year the IF of the journals is calculated as the average of total citations received by papers published in the previous two years. 
#      All metrics are calculated under two scenarios: only main text citations count; and both main text and supplementary citations count.

#The simulation is replicated 10 times.

#set directory with parameter data and functions
files_directory<-'~/Dropbox/CitationFieldEcologists/Codes/POLISHED CODE/'

#load simulation functions
source(paste0(files_directory,'SimulationFunctions.R'))

#Load Data
#references cited per journal in main and sm
refsInfo<-read.csv(paste0(files_directory,'AllRefsInfo.csv'))

#Authors Info #Authors_JournalsProportionsAndPubRate_newPhytEcol
AuthorsINFO<-read.csv(paste0(files_directory, 'Authors_JournalsProportionsAndPubRate.csv'))
#delete authors with zero proportions per category
AuthorsINFO<-AuthorsINFO[apply(AuthorsINFO[,c('Multidisciplinary', 'BiodivCons', 'FieldEcology', 'Ecology', 'Macroecology')], 1, sum)>0,]

citMatrix_m_abs<-read.table(paste0(files_directory,'CrossCitationMatrixMAIN.txt'))
citMatrix_m<-citMatrix_m_abs/rowSums(citMatrix_m_abs)
citMatrix_sm_abs<-read.table(paste0(files_directory,'CrossCitationMatrixSUP.txt'))
citMatrix_sm<-citMatrix_sm_abs/rowSums(citMatrix_sm_abs)

#n papers per journal #NpapersPerJournal_NewPhy_JAnimEcol_Ecology
NpapersPerJournal<-read.csv(paste0(files_directory, 'NpapersPerJournal.csv'))

#proportion of papers cited in a given year
propY_m<-read.csv(paste0(files_directory, 'PropCitPerYearM_loess.csv'))
propY_sm<-read.csv(paste0(files_directory, 'PropCitPerYearSM_loess.csv'))

#Acronmys of journals
Jacr=as.character(NpapersPerJournal[, 'Jacr'])
Jacr_field=as.character(NpapersPerJournal[NpapersPerJournal$Type=='FieldEcology', 'Jacr'])
Jacr_macro=as.character(NpapersPerJournal[NpapersPerJournal$Type=='Macroecology', 'Jacr'])
Jacr_multi=as.character(NpapersPerJournal[NpapersPerJournal$Type=='Multidisciplinary', 'Jacr'])
Jacr_cons=as.character(NpapersPerJournal[NpapersPerJournal$Type=='Conservation', 'Jacr'])
Jacr_ecol=as.character(NpapersPerJournal[NpapersPerJournal$Type=='Ecology', 'Jacr'])

#Run Simulation
S<-RunSimulation()

save(S, file=paste0(files_directory, 'SimulationOutput10reps.Rdata'))

load(paste0(files_directory, 'SimulationOutput10reps.Rdata'))

JournalsIF_m=S[[1]] #Journals' Impact factors accounting for citations in main text only
JournalsIF_sm=S[[2]] #Journals' Impact factors accounting for citations in main text and supplementary materials
AuthorsCit_m=S[[3]] #Authors citations per year accounting for citations in main text only
AuthorsCit_sm=S[[4]] #Authors citations per year accounting for citations in main text and supplementary materials
AuthorsHIndex_m=S[[5]] #Authors H-Index per year accounting for citations in main text only
AuthorsHIndex_sm=S[[6]] #Authors H-Index per year accounting for citations in main text and supplementary materials
PublicationsInfo<-S[[7]] #list with publication info, 1 dataframe per year
AuthorsInfo<-S[[8]] #Authors info

#replicate functions n times
#S=ReplicateSimulation(nRep=10)

#Average results from replicates
#JournalsIF_m=aveReplicates(S, 1, 'Journal') #Journals IF per year main references
#JournalsIF_sm=aveReplicates(S, 2, 'Journal') #Journals IF per year main and SM references
#AuthorsCit_m=aveReplicates(S, 3, 'Author')  
#AuthorsCit_sm=aveReplicates(S, 4, 'Author')
#AuthorsHIndex_m=aveReplicates(S, 5, 'Author')
#AuthorsHIndex_sm=aveReplicates(S, 6, 'Author')
#out<-list()
#for (i in 1:10) {out[[i]]<-do.call('rbind', S[[i]][[7]])}
#PublicationsInfo<-do.call('rbind', out)
#AuthorsInfo<-S[[1]][[8]]

