# Last update: 05.01.2008

### This function calls RAxML (Stamatakis 2006) from within R to estimate topology and branch lengths or do non-parametric bootstrapping. Author: C. Heibl

raxml <- function(seq, runs = 10, bs = FALSE, graph = FALSE, outgroup, path = FALSE, optimize = TRUE, clear = FALSE, rm.ident = FALSE, file = "fromR"){
	
	# fix path:
	if (!is.character(path)){
		path <- system("locate raxmlHPC", intern=TRUE)
		path <- path[1]
		path <- gsub("raxmlHPC", "", path)
	}
	cat(paste("\npath was set to", path, "\n"))	
	
	# definitions:
	R.wd <- getwd()
	setwd(path)
	seq <- delete.empty.cells(seq)
	# remove identical sequences
	if (rm.ident) seq <- identical.seq(seq, delete = TRUE)
	
	# write phylip file
	write.phy(seq, file)
	
	##########################################################
	# 	test initial rearrangement
	##########################################################
	
	if (optimize){
		system(paste("./raxmlHPC -y -s ", file ,".phy -m GTRCAT -n ST0", sep=""))
		system(paste("./raxmlHPC -y -s ", file ,".phy -m GTRCAT -n ST1", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -y -s ", file ,".phy -m GTRCAT -n ST2", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -y -s ", file ,".phy -m GTRCAT -n ST3", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -y -s ", file ,".phy -m GTRCAT -n ST4", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST0 -n FI0", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST1 -n FI1", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST2 -n FI2", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST3 -n FI3", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -i 10 -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST4 -n FI4", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST0 -n AI0", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST0 -n AI1", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST0 -n AI2", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST0 -n AI3", sep=""), show.output.on.console = FALSE)
		system(paste("./raxmlHPC -f d -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST0 -n AI4", sep=""), show.output.on.console = FALSE)
		x <- scan("RAxML_info.FI0", what="c", sep="", quiet = TRUE)
		FI0 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.FI1", what="c", sep="", quiet = TRUE)
		FI1 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.FI2", what="c", sep="", quiet = TRUE)
		FI2 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.FI3", what="c", sep="", quiet = TRUE)
		FI3 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.FI4", what="c", sep="", quiet = TRUE)
		FI4 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.AI0", what="c", sep="", quiet = TRUE)
		AI0 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.AI1", what="c", sep="", quiet = TRUE)
		AI1 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.AI2", what="c", sep="", quiet = TRUE)
		AI2 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.AI3", what="c", sep="", quiet = TRUE)
		AI3 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.AI4", what="c", sep="", quiet = TRUE)
		AI4 <- x[grep("Final", x)-1]
		FIXED <- mean(as.numeric(c(FI0, FI1, FI2, FI3, FI4)))
		AUTO <- mean(as.numeric(c(AI0, AI1, AI2, AI3, AI4)))
		# set i
		if (AUTO > FIXED) {
			i <- x[grep("setting", x) + 1]
			i <- as.numeric(gsub(",", "", i))
		}
		else i <- 10
	
	##########################################################
	# 	number of categories
	##########################################################
		x <- paste("./raxmlHPC -f d -i ", i, " -m GTRMIX -s ", file ,".phy -t RAxML_parsimonyTree.ST", sep="")
		system(paste(x, "0 -c 10 -n C10_0", sep=""), show.output.on.console = FALSE)
		system(paste(x, "1 -c 10 -n C10_1", sep=""), show.output.on.console = FALSE)
		system(paste(x, "2 -c 10 -n C10_2", sep=""), show.output.on.console = FALSE)
		system(paste(x, "3 -c 10 -n C10_3", sep=""), show.output.on.console = FALSE)
		system(paste(x, "4 -c 10 -n C10_4", sep=""), show.output.on.console = FALSE)
		system(paste(x, "0 -c 40 -n C40_0", sep=""), show.output.on.console = FALSE)
		system(paste(x, "1 -c 40 -n C40_1", sep=""), show.output.on.console = FALSE)
		system(paste(x, "2 -c 40 -n C40_2", sep=""), show.output.on.console = FALSE)
		system(paste(x, "3 -c 40 -n C40_3", sep=""), show.output.on.console = FALSE)
		system(paste(x, "4 -c 40 -n C40_4", sep=""), show.output.on.console = FALSE)
		system(paste(x, "0 -c 55 -n C55_0", sep=""), show.output.on.console = FALSE)
		system(paste(x, "1 -c 55 -n C55_1", sep=""), show.output.on.console = FALSE)
		system(paste(x, "2 -c 55 -n C55_2", sep=""), show.output.on.console = FALSE)
		system(paste(x, "3 -c 55 -n C55_3", sep=""), show.output.on.console = FALSE)
		system(paste(x, "4 -c 55 -n C55_4", sep=""), show.output.on.console = FALSE)
		x <- scan("RAxML_info.C10_0", what="c", sep="", 
			quiet = TRUE)
		C10_0 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.C10_1", what="c", sep="", 
			quiet = TRUE)
		C10_1 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.C10_2", what="c", sep="", 
			quiet = TRUE)
		C10_2 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.C10_3", what="c", sep="", 
			quiet = TRUE)
		C10_3 <- x[grep("Final", x)-1]
		x <- scan("RAxML_info.C10_4", what="c", sep="", 
			quiet = TRUE)
		C10_4 <- x[grep("Final", x)-1]
		C10 <- mean(as.numeric(c(C10_0, C10_1, C10_2, C10_3, 			C10_4)))
	
	x <- scan("RAxML_info.C40_0", what="c", sep="", quiet = TRUE)
	C40_0 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C40_1", what="c", sep="", quiet = TRUE)
	C40_1 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C40_2", what="c", sep="", quiet = TRUE)
	C40_2 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C40_3", what="c", sep="", quiet = TRUE)
	C40_3 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C40_4", what="c", sep="", quiet = TRUE)
	C40_4 <- x[grep("Final", x)-1]
	C40 <- mean(as.numeric(c(C40_0, C40_1, C40_2, C40_3, C40_4)))
	
	x <- scan("RAxML_info.C55_0", what="c", sep="", quiet = TRUE)
	C55_0 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C55_1", what="c", sep="", quiet = TRUE)
	C55_1 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C55_2", what="c", sep="", quiet = TRUE)
	C55_2 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C55_3", what="c", sep="", quiet = TRUE)
	C55_3 <- x[grep("Final", x)-1]
	x <- scan("RAxML_info.C55_4", what="c", sep="", quiet = TRUE)
	C55_4 <- x[grep("Final", x)-1]
	C55 <- mean(as.numeric(c(C55_0, C55_1, C55_2, C55_3, C55_4)))
	
	C <- if (FIXED > AUTO) c(C10, FIXED, C40, C55) else c(C10, AUTO, C40, C55)
	catnum <- c(10, 25, 40, 55)
	DF <- cbind(catnum, C)
	DF <- DF[order(DF[,2], decreasing = TRUE),]
	numcat <- DF[1,1]
	
	} # end of OPTIMIZE
	
	##########################################################
	# 	Find BEST TREE	
	##########################################################
	
	if (bs == FALSE){
		if (optimize) call.phy <- paste("./raxmlHPC -f d -i ", i, 
			" -c ", numcat, " -m GTRMIX -s ", file ,".phy -# ", 			runs, " -n ", file, ".tre", sep="")
		else call.phy <- paste("./raxmlHPC -f d -m GTRMIX -s ", 			file , "-p", 666661, ".phy -# ", runs, " -n ", file, ".tre", sep="")
		system(call.phy)
		if (optimize){
			cat(paste("\nE(logLik) with fixed rearrangement settings:", FIXED))
			cat(paste("\nE(logLik) with automatic rearrangement settings:", AUTO))
		}
		cat("\nRAxML was called as follows:")
		cat(call.phy) 
		x <- scan(paste(path, "RAxML_info.", file, ".tre",sep=""), 			what="c", sep="", quiet = TRUE)
		i <-as.integer(gsub(":", "", x[grep("Best", x) + 5]))
		BEST.TREE <- read.tree(paste("RAxML_result.", file, 			".tre.RUN.", i, sep=""))
		setwd(R.wd)
		if (graph) plot(BEST.TREE)
		TREENAME <- paste("tre.", file, sep="")
		write.tree(BEST.TREE, file=TREENAME, multi.line = FALSE)
		cat("\n\nThe best tree obtained was printed as ", TREENAME, 			" to your working directory", sep="")
		if (clear){
			system(paste("mkdir", file))
			PATH <- paste(R.wd, "/", file, "/", sep = "")
			setwd(path)
			system(paste("mv ", file, ".phy ", PATH, file, ".phy", 				sep = ""))
			system(paste("mv RAxML_info.", file, ".tre ", PATH, 				"RAxML_info.", file, ".tre", sep = ""))
			for (i in 0:(runs-1)){
				system(paste("mv RAxML_log.", file, ".tre.RUN.", i, 					" ", PATH, "RAxML_log.", file, ".tre.RUN.", i, 					sep = ""))
				system(paste("mv RAxML_parsimonyTree.", file, 					".tre.RUN.", i, " ", PATH, 						"RAxML_parsimonyTree.", 							file, ".tre.RUN.", i, sep = ""))
				system(paste("mv RAxML_result.", file, ".tre.RUN.", 					i, " ", PATH, "RAxML_result.", file, ".tre.RUN.", 					i, sep = ""))
			}
			cat("\n\nThe entire RAxML output has been stored in ", PATH, "\n\n", 			 sep="")
		}
		else cat("\n\nA log file, parsimony starting trees, and the other resulting trees are stored in ", path, "\n\n",sep="")
		BEST.TREE
	}
		
	##########################################################
	# NON-PARAMETRIC BOOTSTRAPPING
	##########################################################	
	else {
		if (optimize){
			call.phy <- paste("./raxmlHPC -f d -i ", i, " -c ", 				numcat, " -m GTRMIX -s ", file ,".phy -# 10 -n ", 				file, ".tre", sep="") # get best tree
			call.bs <- paste("./raxmlHPC -f d -i ", i, " -c ", 				numcat, " -m GTRMIX -s ",file ,".phy -# ", runs, 
				" -b ", bs," -n ", file, ".boot", sep="")#bootstrap
		}
		else {
			call.phy <- paste("./raxmlHPC -f d -i ", i, " -c ", 				numcat, " -m GTRMIX -s ", file ,".phy -# 10 -n ", 				file, ".tre", sep="") # get best tree
			call.bs <- paste("./raxmlHPC -f d -i ", i, " -c ", 				numcat, " -m GTRMIX -s ",file ,".phy -# ", runs, 
				" -b ", bs," -n ", file, ".boot", sep="")#bootstrap
			}
		system(call.phy)
		system(call.bs)
		cat("\nRAxML was called as follows:")
		cat(paste("\n", call.phy, sep = ""))
		cat(paste("\n", call.bs, sep = ""))
		# 1. Identify best tree
		x <- scan(paste(path, "RAxML_info.", file, ".tre",sep=""), 			what="c", sep="", quiet = TRUE)
		i <-as.integer(gsub(":", "", x[grep("Best", x) + 5]))
		# 2. Read in best tree
		BEST.TREE <- read.tree(paste("RAxML_result.", file, 			".tre.RUN.", i, sep=""))
		# 3. Root best tree
		BEST.TREE <- root(BEST.TREE, outgroup)
		# 4. schreibe ihn in RAxml: rootedbesttreefromR 
		write.tree(BEST.TREE, file="rootedbesttreefromR", 			multi.line=FALSE)
		system(paste("./raxmlHPC -f b -m GTRCAT -s ", file, 
			".phy -z RAxML_bootstrap.", file , 
			".boot -t rootedbesttreefromR -n", file, ".bipart", 			sep=""))
		# 5. Read supported tree
		SUPPORT.TREE <- read.tree(paste("RAxML_bipartitions.", 			file, ".bipart", sep=""))
		setwd(R.wd)
		TREENAME <- paste("bootstrap.", file, sep="")
		write.tree(SUPPORT.TREE, TREENAME)
		if (graph == TRUE){
			plot(SUPPORT.TREE, no.margin=TRUE, y.lim=length				(SUPPORT.TREE$tip.label))
			bootvalues(list(SUPPORT.TREE$node.label), minbs=50)
		}	
	}
	setwd(R.wd)
}








### This functions calculates bipartions using RAxML algorithms and returns them to R. Author: C. Heibl
 
raxml.bipart <- function(seq, tree, boot, file, path="/Applications/RAxML-VI-HPC-2.2.3/")
	{
	if(!all(names(seq) %in% tree$tip.label))
		stop("seq and tree are not identical!")
	RWD <- getwd()
	setwd(path)
	# transfer files		
	write.dna.phylip(seq, "SEQ")
	write.tree(tree, "/Applications/RAxML-VI-HPC-2.2.3/TREE")
	write.multi.tree(boot, "/Applications/RAxML-VI-HPC-2.2.3/BOOT", multi.line=FALSE)
	
	# calculate bipartions
	system("./raxmlHPC -f b -m GTRCAT -s SEQ.phylip -z BOOT -t TREE -n BIPART")
	
	# get bipartitions
	output.tree <- read.tree("RAxML_bipartitions.BIPART")
	
	# delete files in the raxml folder
	system("rm SEQ.phylip")
	system("rm TREE")
	system("rm BOOT")
	system("rm RAxML_info.BIPART")
	system("rm RAxML_bipartitions.BIPART")
	
	setwd(RWD)
	invisible(output.tree)
	
	}
	
	
	
	
	
	




