require(scales)
require(plyr)
require(partitions)

RunSimulation<-function(
					nYears=10
					) 
{

nAuthors=nrow(AuthorsINFO)

#generate journals
Journals=data.frame(Journal=NpapersPerJournal$Jacr, Type=NpapersPerJournal$Type)
Journals$IF<-0
Journals$IF_all<-0

#generate authors
Authors<-AuthorsINFO
Authors$Author<-paste0('A',1:sum(nAuthors))
names(Authors)[which(names(Authors)=='BiodivCons')]<-'Conservation'
Authors$Cit=0
Authors$Cit_all=0

#create outputs to fill
JournalsIF_m=data.frame(Journal=Journals$Journal, Type=Journals$Type)
JournalsIF_sm=data.frame(Journal=Journals$Journal, Type=Journals$Type)
AuthorsCit_m=data.frame(Author=Authors$Author)
AuthorsCit_sm=data.frame(Author=Authors$Author)
AuthorsHIndex_m=data.frame(Author=Authors$Author)
AuthorsHIndex_sm=data.frame(Author=Authors$Author)
Publications=list()#to store publications from previous years

#loop through years
	for (y in 1:nYears) {

		print(paste('year =', y))
		
		#calculate number of papers published this year sampling from a poisson distribution
		Authors$nPapers<-rpois(nrow(Authors), Authors$Rate)

		#update Papers with the new papers published
		Papers=data.frame(Author=rep(Authors$Author, times=Authors$nPapers), TotalCit=0, TotalCit_SM=0)
		Papers$Paper<-paste0('MS', 1:nrow(Papers), 'y', y)
		Papers<-merge(Papers, Authors[,c('Author','Cit', 'Cit_all', 'Multidisciplinary', 'Conservation', 'FieldEcology', 'Ecology', 'Macroecology')], by='Author', all.x=TRUE) #add author cit
		
		#assign journals to the new papers
			for (p in 1:nrow(Papers)) {
			#sample journal category
			Papers[p,'Type']<-sample(c('Multidisciplinary', 'Conservation', 'FieldEcology', 'Ecology', 'Macroecology'), 1, prob=Papers[p, c('Multidisciplinary', 'Conservation', 'FieldEcology', 'Ecology', 'Macroecology')], replace=FALSE)
			#sample journal based on n of papers per journal in that category
			if(Papers[p,'Type']=='FieldEcology') {Papers[p, 'Journal']<-sample(NpapersPerJournal$Jacr[NpapersPerJournal$Jacr %in% Jacr_field], 1, prob=NpapersPerJournal[NpapersPerJournal$Jacr %in% Jacr_field,'Freq'])}
			if(Papers[p,'Type']=='Macroecology') {Papers[p, 'Journal']<-sample(NpapersPerJournal$Jacr[NpapersPerJournal$Jacr %in% Jacr_macro], 1, prob=NpapersPerJournal[NpapersPerJournal$Jacr %in% Jacr_macro,'Freq'])}
			if(Papers[p,'Type']=='Conservation') {Papers[p, 'Journal']<-sample(NpapersPerJournal$Jacr[NpapersPerJournal$Jacr %in% Jacr_cons], 1, prob=NpapersPerJournal[NpapersPerJournal$Jacr %in% Jacr_cons,'Freq'])}
			if(Papers[p,'Type']=='Multidisciplinary') {Papers[p, 'Journal']<-sample(NpapersPerJournal$Jacr[NpapersPerJournal$Jacr %in% Jacr_multi], 1, prob=NpapersPerJournal[NpapersPerJournal$Jacr %in% Jacr_multi,'Freq'])}
			if(Papers[p,'Type']=='Ecology') {Papers[p, 'Journal']<-sample(NpapersPerJournal$Jacr[NpapersPerJournal$Jacr %in% Jacr_ecol], 1, prob=NpapersPerJournal[NpapersPerJournal$Jacr %in% Jacr_ecol,'Freq'])}
		}
		
		Papers<-merge(Papers, Journals[, c('Journal', 'IF', 'IF_all')], by='Journal')	

		#assign number of references in main and supplementary materials by sampling from the data collected
		#these are the citations that every paper makes
		for (i in 1:nrow(Papers)) {
			Papers$J_Cits[i]<-sample(refsInfo[refsInfo$Jacr==as.character(Papers$Journal[i]),'nRefs'], 1, replace=TRUE)
			Papers$J_SM_Cits[i]<-sample(refsInfo[refsInfo$Jacr==as.character(Papers$Journal[i]),'nSupRefs'], 1, replace=TRUE)
		}	
		
		#calculate the total number of citations based on the total number of references
		TotalCit<-tapply(Papers$J_Cits, Papers$Journal, sum)
		TotalCitSM<-tapply(Papers$J_SM_Cits, Papers$Journal, sum)

		Papers$PapersCit<-0
		Papers$PapersCit_SM<-0
		Publications[[y]]<-Papers	

			#Calculate new citations for papers published in previous years
			for (i in 1:y) { #For all years before present, add new citations to published papers and re-calculate total citations
				
				#calculate the number of times the previous n of citations should be multipled
				pYm=propY_m[propY_m$Year==2017-(y-i), 'Prop']
				pYsm=propY_sm[propY_sm$Year==2017-(y-i), 'Prop']	

				#these are the citations that each journal makes to papers in this year
				TotalCit2<-TotalCit*pYm
				TotalCitSM2<-TotalCitSM*pYsm			

				#Calculate the proportion of papers published in the simulation by journal compared to the total in reality
				tab<-table(Papers$Journal)
				PapersPublished=data.frame(Jacr=names(tab), N=as.numeric(table(Papers$Journal)))
				PapersPublished<-merge(PapersPublished, NpapersPerJournal[,c('Jacr', 'Freq')], by='Jacr', all.x=TRUE)
				PapersPublished$Prop<-PapersPublished$N/PapersPublished$Freq
				
				#Adjust considering citations coming from papers outside the simulation
				PapersPublished2<-merge(PapersPublished, data.frame(Jacr=names(TotalCit2), Cit=as.numeric(TotalCit2), CitSM=TotalCitSM2), by='Jacr', all.x=TRUE)
				PapersPublished2$Cit<-PapersPublished2$Cit / PapersPublished2$Prop
				PapersPublished2$CitSM<-PapersPublished2$CitSM / PapersPublished2$Prop

				TotalCit2<-PapersPublished2$Cit; names(TotalCit2)<-PapersPublished2$Jacr
				TotalCitSM2<-PapersPublished2$CitSM; names(TotalCitSM2)<-PapersPublished2$Jacr

				for (jg in 1:length(Jacr)) { #loop through journals (for each journal calculates the amount of citations made)

					Jgiver=as.character(Jacr[jg])
					
					Cit_main <- TotalCit2[Jgiver]
					Cit_sm <- TotalCitSM2[Jgiver]

					if(is.na(Cit_main)) {next} #if nobody published in this journal this year, skip to the next

					#calculate how the citations made are distributed to the other journals using the estimated proportions
					Cit_main<-round(citMatrix_m[Jgiver,] * Cit_main)
					names(Cit_main)<-colnames(citMatrix_m)
					Cit_sm<-round(citMatrix_sm[Jgiver,] * Cit_sm)
					names(Cit_sm)<-colnames(citMatrix_sm)

					for (jr in 1:length(Jacr)) { #loop through all papers of each journals (Jreceiver) and distribute the citations produced by journal (J)
						Jreceiver=Jacr[jr] #journal that receives citations
						citDistribute_main=Cit_main[Jreceiver][1,1]
						citDistribute_sm=Cit_sm[Jreceiver][1,1]

						#update citations considering that many citations go to papers that are not in the simulations
						#multiply citations from X to Y by the proportion of papers in Y over the total papers published in Y in reality
						citDistribute_main = round(citDistribute_main * PapersPublished[PapersPublished$Jacr==Jreceiver,'Prop'])
						citDistribute_sm = round(citDistribute_sm * PapersPublished[PapersPublished$Jacr==Jreceiver,'Prop'])

						#extract papers in the journal that receives the citations
						papersJreceiver<-Publications[[i]][Publications[[i]]$Journal==Jreceiver, ]

						if (i == y) {#if adding citations to papers published this year
						
								if(citDistribute_main>0) {
								#if cit > t, parts matrix becomes huge it gets stuck. 
								#So I do it for a smaller number and multiplicate the matrix for the same. Once rounded the sum of citations per column isn't exactly the same original number but errors are limited. However it increase the number of zeros
								t=20
								if(citDistribute_main<t) {pp_m<-parts(citDistribute_main)} else {pp_m<-round(parts(round(citDistribute_main*(t/citDistribute_main)))*(citDistribute_main/t))}
 								pp_m<-pp_m[,apply(pp_m, 2, function(x){sum(x!=0)})<=nrow(papersJreceiver)] #only retain combinations where clump of citations to distribute is < papers (pp_m <= number of papers)
 								if(class(pp_m)=='matrix') {pp_m<-pp_m[,sample(ncol(pp_m), 1)]}#pick a random column (if it is a matrix)
 								pp_m<-pp_m[pp_m!=0] #removes the zeros
 								if(length(pp_m)>1) {pp_m<-pp_m[sample(length(pp_m), length(pp_m), replace=FALSE)]} #reshuffle if >1
								#assign the citations to random papers in the journals
 								sr_m=sample(nrow(papersJreceiver), length(pp_m))#sample rows for main text citations
 								papersJreceiver[sr_m,'PapersCit']<-pp_m + papersJreceiver[sr_m,'PapersCit'] #sum to those already received
								}
						
								if(citDistribute_sm>0) {
								#if cit > t, parts matrix becomes huge it gets stuck. So I do it for a smaller number and multiplicate the matrix for the same. Once rounded the sum of citations per column isn't exactly the same original number but errors are limited.  However it increase the number of zeros
								if(citDistribute_sm<t) {pp_sm<-parts(citDistribute_sm)} else {pp_sm<-round(parts(round(citDistribute_sm*(t/citDistribute_sm)))*(citDistribute_sm/t))}
 								pp_sm<-pp_sm[,apply(pp_sm, 2, function(x){sum(x!=0)})<=nrow(papersJreceiver)] #only retain combinations where clump of citations to distribute is < papers (pp_m <= number of papers)
 								if(class(pp_sm)=='matrix') {pp_sm<-pp_sm[,sample(ncol(pp_sm), 1)]}#pick a random column (if it is a matrix)
 								pp_sm<-pp_sm[pp_sm!=0] #removes the zeros
 								if(length(pp_sm)>1) {pp_sm<-pp_sm[sample(length(pp_sm), length(pp_sm), replace=FALSE)]} #reshuffle if >1
								#assign the citations to random papers in the journals
 								sr_sm=sample(nrow(papersJreceiver), length(pp_sm))#sample rows for sm citations
 								papersJreceiver[sr_sm,'PapersCit_SM']<-pp_sm + papersJreceiver[sr_sm,'PapersCit_SM'] #sum to those already received
 								}

 							} else { #if papers published in previous years
							#sample from a poisson distribution centered to the averate citations received in the previous year
							CitPois_m<-rpois(nrow(papersJreceiver), (papersJreceiver$TotalCit/i)+1) #+1 so papers that were never cited have a chance to be cited
							CitPois_m<-round(citDistribute_main * (CitPois_m / sum(CitPois_m))) #normalize to make the sum equal to the number of citations to distribute
							papersJreceiver$PapersCit<-CitPois_m+papersJreceiver$PapersCit
						
							CitPoisSM<-rpois(nrow(papersJreceiver), (papersJreceiver$TotalCit_SM/i)+1) #+1 so papers that were never cited have a chance to be cited
							CitPoisSM<-round(citDistribute_sm * (CitPoisSM / sum(CitPoisSM))) #normalize to make the sum equal to the number of citations to distribute
							papersJreceiver$PapersCit_SM<-CitPoisSM+papersJreceiver$PapersCit_SM
 							}

						Publications[[i]][Publications[[i]]$Journal==Jreceiver, ]<-papersJreceiver #reinsert the journal in the list
						} #closes loop through papers (Jreceiver)

					} #closes loop through papers (Jgiver)

				#update total cit per publications with the new citations
				Publications[[i]]$TotalCit<-Publications[[i]]$TotalCit + Publications[[i]]$PapersCit
				Publications[[i]]$TotalCit_SM<-Publications[[i]]$TotalCit_SM + Publications[[i]]$PapersCit_SM

				}	#closes citation update (i)
	

			#Update authors' citations
			#calculate the total number of citations authors have received		
			t<-lapply(Publications, tapp) #sum of the citations per author dataframe in the list
			AuthorCit<-tapply(unlist(t), names(unlist(t)), sum) #sum total
			
			AuthorCit=data.frame(Author=names(AuthorCit), newCit=AuthorCit)	
			AuthorCit$newCit[is.na(AuthorCit$newCit)]<-0

			t<-lapply(Publications, tapp2)
			AuthorCit_SM<-tapply(unlist(t), names(unlist(t)), sum)

			AuthorCit_SM=data.frame(Author=names(AuthorCit_SM), newCitSM=AuthorCit_SM)	
			AuthorCit_SM$newCitSM[is.na(AuthorCit_SM$newCitSM)]<-0

			Authors<-merge(Authors, AuthorCit, by='Author', all.x=TRUE)
			Authors<-merge(Authors, AuthorCit_SM, by='Author', all.x=TRUE)
			Authors$newCit[is.na(Authors$newCit)]<-0
			Authors$newCitSM[is.na(Authors$newCitSM)]<-0
			Authors$Cit<-Authors$Cit+Authors$newCit #update cit main text
			Authors$Cit_all<-Authors$Cit_all + Authors$newCit + Authors$newCitSM #update all cit

			Authors<-Authors[,-which(names(Authors) %in% c('newCit', 'newCitSM'))]	

			#Calculate H-Index
			auths<-unique(Authors$Author)
			outH<-data.frame(Author=auths, HIndex=NA, HIndex_all=NA)
			for (i in 1:nrow(outH)) {
				auth<-auths[i]
				l<-lapply(Publications, extractAuthor, author=auth)
				df <- ldply(l, data.frame)
				if(nrow(df)>0) { #if there is at least one publication
				#main
				df2<-aggregate(Cit ~ Paper, df, sum)
				df2$rank<-1:nrow(df2)
				HIndex<-nrow(df2[df2$Cit>=df2$rank,])
				#main + supplementary materials
				df3<-aggregate(Cit_all ~ Paper, df, sum)
				df3$rank<-1:nrow(df3)
				HIndex2<-nrow(df3[df3$Cit_all>=df3$rank,])
				} else {HIndex<-0; HIndex2<-0} #if the author has no publications
				outH[outH$Author==auth, 'HIndex']<-HIndex
				outH[outH$Author==auth, 'HIndex_all']<-HIndex2
			}
			
			#Recalculate journals' Impact Factor every 2 years
			if (y>=3) {#from the third year
		
				PubIF<-rbind(Publications[[y-2]], Publications[[y-1]])

				PubIF$TotalCit_SM[is.na(PubIF$TotalCit_SM)]<-0

				aveCitJ<-tapply(PubIF$TotalCit, PubIF$Journal, mean)
				aveCitJ_all<-tapply(PubIF$TotalCit + PubIF$TotalCit_SM, PubIF$Journal, mean)	

				newIF=data.frame(Journal=names(aveCitJ), IF=aveCitJ, IF_all=aveCitJ_all)
				Journals<-Journals[,-grep('IF', names(Journals))]
				Journals<-merge(Journals, newIF, by='Journal', all.x=TRUE)
				Journals$IF<-ifelse(is.na(Journals$IF), 0, Journals$IF)
				Journals$IF_all<-ifelse(is.na(Journals$IF_all), 0, Journals$IF_all)
				} #closes if IF y2	
	
	

	#store info
	JournalsIF_m<-merge(JournalsIF_m, Journals[,c('Journal', 'IF')], by='Journal', all.x=TRUE) #add IF
	names(JournalsIF_m)[3:ncol(JournalsIF_m)]<-paste0(rep('y'),1:(ncol(JournalsIF_m)-2))
	JournalsIF_sm<-merge(JournalsIF_sm, Journals[,c('Journal', 'IF_all')], by='Journal', all.x=TRUE) #add IF
	names(JournalsIF_sm)[3:ncol(JournalsIF_sm)]<-paste0(rep('y'),1:(ncol(JournalsIF_sm)-2))	

	AuthorsCit_m<-merge(AuthorsCit_m, Authors[,c('Author', 'Cit')], by='Author', all.x=TRUE) #add total citations
	names(AuthorsCit_m)[2:ncol(AuthorsCit_m)]<-paste0(rep('y'),1:(ncol(AuthorsCit_m)-1))
	AuthorsCit_sm<-merge(AuthorsCit_sm, Authors[,c('Author', 'Cit_all')], by='Author', all.x=TRUE) #add total citations
	names(AuthorsCit_sm)[2:ncol(AuthorsCit_sm)]<-paste0(rep('y'),1:(ncol(AuthorsCit_sm)-1))	

	AuthorsHIndex_m<-merge(AuthorsHIndex_m, outH[,c('Author', 'HIndex')], by='Author', all.x=TRUE) #add total citations
	names(AuthorsHIndex_m)[2:ncol(AuthorsHIndex_m)]<-paste0(rep('y'),1:(ncol(AuthorsHIndex_m)-1))
	AuthorsHIndex_sm<-merge(AuthorsHIndex_sm, outH[,c('Author', 'HIndex_all')], by='Author', all.x=TRUE) #add total citations
	names(AuthorsHIndex_sm)[2:ncol(AuthorsHIndex_sm)]<-paste0(rep('y'),1:(ncol(AuthorsHIndex_sm)-1))	

	} #end year loop

	#clean output for authors
	Authors<-Authors[,c('Author','Rate','Multidisciplinary', 'Conservation', 'FieldEcology', 'Ecology', 'Macroecology')]
	for (i in 1:y) {Publications[[i]]<-Publications[[i]][,c('Author', 'Paper', 'TotalCit','TotalCit_SM','Journal','Type')]}

return(list(JournalsIF_m, JournalsIF_sm, AuthorsCit_m, AuthorsCit_sm, AuthorsHIndex_m, AuthorsHIndex_sm, Publications, Authors))

} #close function

extractAuthor<-function(x, author){x<-x[x$Author==author,]; return(x)}

tapp<-function(X){
	t<-tapply(X$PapersCit, X$Author, sum, na.rm=TRUE); t[is.na(t)]<-0; return(t)
}

tapp2<-function(X){
	t<-tapply(X$PapersCit_SM, X$Author, sum, na.rm=TRUE); t[is.na(t)]<-0; return(t)
}

ReplicateSimulation=function(X, nRep=10) {
	out<-list()
	for (i in 1:10) {
	print(paste('RUN REPLICATE', i))
	out[[i]]<-RunSimulation()
	}
	return(out)
}

aveReplicates<-function(S, outputN=1, colRetain='Journal') {
  tmp1<-as.data.frame(S[[1]][[outputN]])[,paste0('y', 1:10)]
  colRet<-as.data.frame(S[[1]][[outputN]])[,colRetain]
  for (i in 1:10) {
    tmp2<-as.data.frame(S[[i]][[outputN]])[,paste0('y', 1:10)]
    tmp1<-(tmp1+tmp2)
  }
  ave<-tmp1/10
  Ave=cbind(colRet, ave)
  return(Ave)
}