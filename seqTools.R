require("Biostrings");
require("BSgenome");
require("rtracklayer");
require("Rsamtools");
require("ShortRead");
require("chipseq")
require("data.table")
require("parallel")

# global constants
gSigDigits = 3;
gNormReads = 1e6

gGenome = "BSgenome.Mmusculus.UCSC.mm9"

#### plotting constants
gMainSize = 1.5;
gLabSize = 1;
gAxisSize = 0.8;
gPch = 19;
gCex = 0.25
gHeight = 3;
gWidth = 3;

gSdLwd = 0.5;
gSdCol = rgb(0.5, 0.5, 0.5);

gRegLwd = 3;
gRegCol = rgb(0, 0, 0.5);
gRegCex = 2;

gMai = c(0.5, 0.5, 0.5, 0.25);
gMgp = c(1.5, 0.5, 0); 

memLog =  function(header = c(NA), file = c(paste0("memory.log.", make.names(Sys.time()), ".txt")), append = c(T), top = c(10), ...) {

	#Debug
	#file = "Bistools.log.debug.txt"; append = F; top = 10; header = NA

	#write header to file or file name
	cat("\n", file = file, append = append)
	if (!is.na(header)) write(header, file = file, append = append) else write(file, file = file, append = append)
	g = gc()
	suppressWarnings(write.table(g, file = file, append = T))
	#m = memory.profile()
	#suppressWarnings(write.table(m, file = file, append = T))

	#sort largest objects in memory by size
	s = sort(sapply(ls(envir = .GlobalEnv),function(x){object.size(get(x))}), decreasing = T) 

	if (!is.na(top) & top > length(s)) top = length(s)
	s = s[1:top]

	write(paste("Top", top, " memory objects"), file = file, append = T)
	suppressWarnings(write.table(s, file = file, append = T, col.names = F))
}

#annotates a bisulfite coverage file
annotBisCov = function(covFile = c(NA), cov = c(NA), txDb, prefix, outFile = if (!is.blank(covFile)) gsub(".cov$|.cov.gz$", ".annot.cov.gz", basename(covFile)) else "annot.cov.gz", removeChrs = c(NA), org = c(NA), keytype = c("ENTREZID")) {

	#Debug	covFile = paste0("zcat < ", sampleCovFile, ".gz"); txDb = tx; prefix = "GRCh37"; outFile = sampleAnnotFile; removeChrs = NA;

	if (is.na(cov) & !is.na(covFile)) cov = fread(covFile, header = T) else stop("No coverage cov or covFile")
	setnames(cov, "rname", "chr")
	setkey(cov, chr, pos)

	if (!is.na(removeChrs)) cov = subset(cov, !grepl(removeChrs, chr))

	covR = GRanges(seqnames = cov$chr, ranges = IRanges(start = cov$pos, end = cov$pos))
	annot = annotBedtoTxDbGene(bed = covR, tx = txDb, prefix = prefix)

	#convert annotation to data table
	annot = as.data.table(as.data.frame(annot))
	setnames(annot, c("seqnames", "start"), c("chr", "pos"))
	setkey(annot, chr, pos)
	gc()

	if (!is.na(org)) {	
		symbol = rep(NA, dim(annot)[1])
		symbol[!is.na(as.character(annot[[paste0(prefix, ".tsKg")]]))] = select(org, keytype = keytype, keys = as.character(annot[[paste0(prefix, ".tsKg")]])[!is.na(as.character(annot[[paste0(prefix, ".tsKg")]]))], columns = c("SYMBOL"))$SYMBOL
		annot[, eval(parse(text = paste0(prefix, ".tsSymbol:=symbol")))]

		symbol = rep(NA, dim(annot)[1])
		symbol[!is.na(as.character(annot[[paste0(prefix, ".tsKg")]]))] = select(org, keytype = keytype, keys = as.character(annot[[paste0(prefix, ".tssKg")]])[!is.na(as.character(annot[[paste0(prefix, ".tssKg")]]))], columns = c("SYMBOL"))$SYMBOL
		annot[, eval(parse(text = paste0(prefix, ".tssSymbol:=symbol")))]
	}

	cov = merge(annot, cov)
	gc()
	
	writeBistoolsCov(cov, outFile)

	invisible(cov)

}

#Compiles sample coverage summarizing by group and enforcing a minimium coverage (minCov) at the group level
compileBisGrpCov = function(samples, grps, covFiles, covFile = c(NA), minCov = c(10), reportSamples = c(T), verbose = c(F), ...) {
	#Debug samples = files$sample; grps = files$sample_group; covFiles = paste0(files$dir, files$covFile); covFile = sampleCovFile; minCov = minCov; reportSamples = F; verbose = T

	bedExt = ".bed";

	if (length(samples) == length(covFiles) & length(samples) == length(grps)) {
		for (i in 1:length(samples)) {

			if (verbose) print(paste("Reading", samples[i], "-", Sys.time()))
			cov = readBistoolsCov(covFiles[i])
	
			#set appropriate names
			cov = setnames(cov, names(cov)[!names(cov) %in% key(cov)], paste(samples[i], names(cov)[!names(cov) %in% key(cov)], sep = "."))
	
			#aggregate
			if (i == 1) covs = cov else covs = merge(covs, cov, all = T)
			rm(cov); gc()
		}

		#Summarize each group
		for (grp in unique(grps)) {
		
			if (verbose) print(paste("Summarizing", grp, "-", Sys.time()))

			grpSamples = samples[grps == grp]
			covs[, eval(parse(text = paste0(grp, ".u:=rowSums(covs[, list(", paste0(grpSamples, ".u", collapse = ","), ")], na.rm = T)"))), ]
			covs[, eval(parse(text = paste0(grp, ".m:=rowSums(covs[, list(", paste0(grpSamples, ".m", collapse = ","), ")], na.rm = T)"))), ]
			covs[, eval(parse(text = paste0(grp, ".meth:=", grp, ".m/(", grp, ".u+", grp, ".m)"))), ]
		}

		if (verbose) print(paste("Filtering for minimum coverage", minCov, "-", Sys.time()))

		#Filter summarized groups
		grpCrit	= list(); #character logical expressions to be evaluated 
		for (grp in unique(grps)) grpCrit[[grp]] = paste0(paste(paste0(grp, c(".u", ".m")), collapse = "+"), ">=", minCov)
		covs = covs[eval(parse(text = paste(grpCrit, collapse = "&")))]

		if (!reportSamples) {
			sampleDataCols = c(".u", ".m", ".meth")
			sampleCols = NULL
			for (col in sampleDataCols) sampleCols = c(sampleCols, paste0(samples, col))
			covs = subset(covs, , !(names(covs) %in% sampleCols))
		}

		#make grp bedGraphs
		for (grp in unique(grps)) writeBistoolsCovToBw(covs, file = paste0(dirname(covFile), "/", grp, ".minCov.", minCov, ".bw"), chrCol = "rname", scoreCol = paste0(grp, ".meth"), ...) 

		if (!is.na(covFile)) {
			writeBistoolsCov(covs, covFile);
			writeBistoolsCovToBed(covs, paste0(covFile, bedExt))
		}
		invisible(covs)

	} else warning("sample length is not the same as coverage files length");	
}

#Merges coverage from several bistools .cov files and imposes a min. cov
#Saves coverage by chromosome and merges chromsomes separately to save memory on merge 
compileBisSampleCovByChr = function(samples, covFiles, covFile = c(paste0("Sample.minCov", minCov, ".minSampleCov", minSampleCov, ".csv")), minCov = c(5), minSampleCov = c(0.8), makeBw = c(T), removeChrs = c(NA), verbose = c(T), ...) {
	#Debug
	#samples = files$sample; covFiles = paste0(files$dir, files$covFile); covFile = sampleCovFile; minCov = minCov; minSampleCov = minSampleCov; makeBw = F; removeChrs = removeChrs; verbose = T; genome = genome; 

	bedExt = ".bed";
	bwExt = ".bw";
	covExt = ".cov"
	gzExt = ".gz"
	
	covs = NULL

	if (length(samples) == length(covFiles)) {
		for (i in 1:length(samples)) {

			cov = readBistoolsCov(covFiles[i])
			
			if(verbose) {
				print(paste(covFiles[i], "-", dim(cov)[1], "lines", "-", Sys.time()))
				print(paste(covFiles[i], "-", minCov, "x coverage at ", sum(cov$u + cov$m >= minCov), " loci"))
			}
				
			#set appropriate names & limit to minCov 
			cov = setnames(cov, names(cov)[!names(cov) %in% key(cov)], paste(samples[i], names(cov)[!names(cov) %in% key(cov)], sep = "."))
	
			#load data into a list of lists of data.tables segregated by chromosome and sample
			chrs = unique(cov$rname)
			for (chr in chrs) covs[[chr]][[samples[i]]] = subset(cov, rname == chr)

			rm(cov);
			if (verbose) {
				print(paste(i, "cov size", format(object.size(covs), units = "Gb"), Sys.time()))
				print(gc())
			} else gc()
		}

		#loop through each chromosome and merge all samples
		for (chr in names(covs)) {
			if (verbose) print(paste("Merging", chr, Sys.time()))
			covs[[chr]] = Reduce(function (x, y) merge(x, y, all = T), covs[[chr]])				
			gc()
		}

		#combine chr data
		covs = rbindlist(covs, fill = T)
		gc()
		print(paste("compiled coverage", Sys.time()))

		#subset for min sample coverage
		covx = subset(covs, , paste0(files$sample, ".u")) + subset(covs, , paste0(files$sample, ".m"))
	
		rowCov = rowSums(covx >= minCov, na.rm = T)
		covs = subset(covs, rowCov >= minSampleCov * length(samples))

		rm(covx); gc()

		print(paste("filtered rows to", minCov, "x coverage at ", sum(rowCov >= minSampleCov * length(samples)), "rows", Sys.time()))

		if (!is.na(covFile)) {
			writeBistoolsCov(covs, covFile);
			writeBistoolsCovToBed(covs, gsub(paste0(covExt, "(", gzExt, ")?$"), bedExt, covFile))
			if (makeBw) for (i in 1:length(samples)) writeBistoolsCovToBw(covs, file = paste0(dirname(covFile), "/", samples[i], ".minCov", minCov, ".minSamples", minSampleCov, bwExt), chrCol = "rname", scoreCol = paste0(samples[i], ".meth"), removeChrs = removeChrs, ...)
		}

	} else warning("sample length is not the same as coverage files length");	

	invisible(covs)
}

#Function to replace NAs in a data table
replaceDtNa = function(dt, value = 0, cols = c(NA)) {

	if (all(is.na(cols))) cols = names(dt)
	for (col in cols) eval(parse(text=paste0("dt[is.na(", col, "), ", col, ":=", value, "]")))
}

#Merges coverage from several bistools .cov files and imposes a min. cov
compileBisSampleCov = function(samples, covFiles, covFile = c(paste0("Sample.minCov", minCov, ".minSampleCov", minSampleCov, ".csv")), minCov = c(5), minSampleCov = c(0.9), makeBw = c(T), removeChrs = c(NA), verbose = c(T), ...) {
	#Debug
	#samples = files$sample; covFiles = paste0(files$dir, files$covFile, ".gz"); covFile = sampleCovFile; minCov = minCov; genome = genome; verbose = T; minSampleCov = minSampleCov; removeChrs = removeChrs;

	bedExt = ".bed";
	bwExt = ".bw";
	covExt = ".cov";
	gzExt = ".gz";
	methSuf = ".meth"
	
	if (length(samples) == length(covFiles)) {
		for (i in 1:length(samples)) {

			cov = readBistoolsCov(covFiles[i])
			
			if(verbose) {
				print(paste(covFiles[i], "-", dim(cov)[1], "lines", "-", Sys.time()))
				print(paste(covFiles[i], "-", minCov, "x coverage at ", sum(cov$u + cov$m >= minCov), " loci"))
			}
	
			#set appropriate names & limit to minCov 
			cov = setnames(cov, names(cov)[!names(cov) %in% key(cov)], paste(samples[i], names(cov)[!names(cov) %in% key(cov)], sep = "."))
	
			#aggregate
			if (i == 1) covs = cov else covs = merge(covs, cov, all = T)

			rm(cov);
			if (verbose) {
				print(paste(i, "cov size", format(object.size(covs), units = "Gb"), Sys.time()))
				print(gc())
			} else gc()
		}
		
			#replace NA values with 0
		replaceDtNa(covs, value = 0, cols = names(covs)[!grepl(methSuf, names(covs))])

		print(paste("compiled coverage", Sys.time()))

		#subset for min sample coverage
		covx = subset(covs, , paste0(files$sample, ".u")) + subset(covs, , paste0(files$sample, ".m"))
	
		rowCov = rowSums(covx >= minCov, na.rm = T)
		covs = subset(covs, rowCov >= minSampleCov * length(samples))

		rm(covx); gc()

		print(paste("filtered rows to", minCov, "x coverage at ", sum(rowCov >= minSampleCov * length(samples)), "rows", Sys.time()))

		if (!is.na(covFile)) {
			writeBistoolsCov(covs, covFile);
			writeBistoolsCovToBed(covs, gsub(paste0(covExt, "(", gzExt, ")?$"), bedExt, covFile))
			if (makeBw) for (i in 1:length(samples)) writeBistoolsCovToBw(covs, file = paste0(dirname(covFile), "/", samples[i], ".minCov", minCov, ".minSamples", minSampleCov, bwExt), chrCol = "rname", scoreCol = paste0(samples[i], ".meth"), removeChrs = removeChrs, ...)

		}

	} else warning("sample length is not the same as coverage files length");	

	invisible(covs)
}

#Function to test if variable exists
is.blank = function(x) return(is.null(x) || length(x) == 0 || all(is.na(x)) || all(x==""))

#Extract sra file
sraExtract = function(sraFile) {

	print(paste("Extracting sra files", Sys.time()))

	sraCmd = "fastq-dump -v --split-3 "
	fq.ext = ".fastq"

	sraCmdFile = paste0(sraCmd, files$dir[i], files$sraFile[i], " -O ", files$dir[i]);
	system(sraCmdFile)

	#bashFile = "sraExt.bash"
	#write(sraCmdFile, file = bashFile, append = if (i == 1) F else T)
	# ssh in and run 'nohup ./sraExt.bash & exit'
}

#Get fastq file names from extracted sra files
getSraFastqNames = function(sraFile) {

	mate1Tag = "_1";
	mate2Tag = "_2";

	dir = dirname(sraFile);
	sraFile = basename(sraFile);

	print(paste("Getting fastq file names", sraFile, Sys.time()))

	fqExt = ".fastq";
	sraExt = ".sra";

	fqFiles = list.files(dir)
	fqMate1 = paste(fqFiles[grepl(paste0(mate1Tag, fqExt, "$"), fqFiles) & grepl(paste0("^", gsub(sraExt, "", sraFile)), fqFiles)], collapse = "|")
	fqMate2 = paste(fqFiles[grepl(paste0(mate2Tag, fqExt, "$"), fqFiles) & grepl(paste0("^", gsub(sraExt, "", sraFile)), fqFiles)], collapse = "|")
	fqFile = paste(fqFiles[grepl(paste0(fqExt, "$"), fqFiles) & grepl(paste0("^", gsub(sraExt, "", sraFile)), fqFiles) & !grepl(paste0(mate1Tag, fqExt, "$"), fqFiles) & !grepl(paste0(mate2Tag, fqExt, "$"), fqFiles)], collapse = "|")

	return(c(fqFile, fqMate1, fqMate2))

}

#call rsync to move files
rsync = function(oldDir, newDir, options = c("-avu")) {
	
	rsyncCall = paste("rsync", options, oldDir, newDir)
	if (system(rsyncCall) == 0) invisible(newDir)
}


#create dir and move files
#oldFiles: pipe (or split) delimited list of files to move
#dir: directory where oldFiles reside
#mvDir: directory to move files to
mvFiles = function(oldFiles, dir, mvDir, split = c("\\|")) {
	#Debug	oldFiles = as.character(files$sraFile[i]); dir = downloadDir; mvDir = paste0(dataDir, files$sample[i], "/"); split = "\\|"
		
	if (!is.blank(oldFiles)) {	

		oldFiles = paste0(dir, strsplit(oldFiles, split)[[1]])
		if (!file.exists(mvDir)) dir.create(mvDir); #create new directory if necessary
		#tryCatch(suppressWarnings(file.rename(oldFiles, paste0(mvDir, basename(oldFiles)))), error = copyDelete(oldFiles, paste0(mvDir, basename(oldFiles))), warning = message("rename not possible"), finally = message("copied and deleted"))
		file.rename(oldFiles, paste0(mvDir, basename(oldFiles)))
		invisible(mvDir)
	}
}

#cp Files dir and move files
cpFiles = function(oldFiles, dir, cpDir, split = c("\\|"), ...) {
	#Debug
	#oldFiles = paste0(files$fqMate1[i], "|", files$fqMate2[i]); dir = files$dir[i]; cpDir = paste0(stageDir, files$sample[i], "/"); split = "\\|";
		
	if (!is.blank(oldFiles)) {	

		oldFiles = paste0(dir, strsplit(oldFiles, split)[[1]])
		if (!file.exists(cpDir)) dir.create(cpDir); #create new directory if necessary
		#tryCatch(suppressWarnings(file.rename(oldFiles, paste0(mvDir, basename(oldFiles)))), error = copyDelete(oldFiles, paste0(mvDir, basename(oldFiles))), warning = message("rename not possible"), finally = message("copied and deleted"))
		file.copy(oldFiles, paste0(cpDir, basename(oldFiles)), ...)
		invisible(cpDir)
	}
}

#copy and delete files, used in case rename.file throws a warning (occurs across file systems) 
copyDelete = function(oldFiles, newFiles) if (file.copy(from = oldFiles, to = newFiles) & file.remove(oldFiles)) return(NA) else return(error("copy and delete failed"))

catFiles = function(files, outFile, dir = c(NA), split = c("\\|"), del = c(F), compress = c(NA), overwrite = c(F)) {

	#Debug	
	#files = paste(paste0(sampleFiles$dir, sampleFiles$collapseCpgFile), collapse = "|"); outFile = paste0(sample, callCollapseExt); dir = NA; split = "\\|"; del = F
	#outFile = paste0(files$sample[i], ".fq", if (gz) ".gz"); files = fqTemp; del = T; split = "\\|";

	#ensure there are files so no error is thrown
	if (!is.blank(files)) {	

		#handle gz compression
		if ((is.na(compress) & grepl("\\.gz$", outFile)) | (!is.na(compress) & compress)) {
			compress = T
			outFile = gsub("\\.gz$", "", outFile)
		} else compress = F

		print(paste("concatenate files", Sys.time()))
	
		if (!is.na(split)) files = unlist(strsplit(files, split))
		if (!is.na(dir)) files = paste0(dir, files) else dir = paste0(dirname(files[1]), "/")

		system(paste("cat", paste(files, collapse = " "), ">", outFile))

		#handle gz compression
		if (compress) {
			system(paste("pigz", if (overwrite) "-f", outFile))
			outFile = paste0(outFile, ".gz")
		}
		
		if (del) file.remove(files)

		invisible(as.character(outFile))
	}
}

formatFastq = function(fqFiles, dir, sampleName = c(NA), outDir = c(NA), split = c("\\|"), gzExt = c(".gz"), fqExt = c(".fastq")) {

	#Debug	fqFiles = files$fqMate1[i]; dir = files$dir[i]; sampleName = paste0(files$sample[i], "_1"); split = "\\|"; gzExt = ".gz"; fqExt = ".fastq";

	#ensure there are files so no error is thrown
	if (!is.null(fqFiles) && !is.na(fqFiles) && fqFiles != "") {	

		print(paste("Formatting fastq files", Sys.time()))
	
		fqFiles = paste0(dir, strsplit(fqFiles, split)[[1]])

		#set outDir
		if (is.na(outDir)) outDir = dir;
		if (!file.exists(outDir)) dir.create(outDir)

		#make output fqFile name
		if (is.na(sampleName)) fqFile = paste0(outDir, gsub(gzExt, "", fqFiles[1])) else fqFile = paste0(outDir, sampleName, fqExt); #combined fastq File name
	
		for (i in 1:length(fqFiles)) {

			#if the file is gzipped
			if (grepl(paste0(gzExt, "$"), fqFiles[i])) {
				system(paste("pigz -cd", fqFiles[i], if (i == 1) ">" else ">>", fqFile)); #decompress files
			} else {
				system(paste("cat", paste0(dir, fqFiles[i]), if (i == 1) ">" else ">>", fqFile))
			}
		}

		return(as.character(basename(fqFile)))
	}
}

#check to see if file is gzipped and if not compress and return name
gzipFile = function(file,  overwrite = c(F), rmIfExists = c(F)) {

	gz.ext = ".gz";

	if (!grepl(paste0(gz.ext, "$"), file) & (!file.exists(paste0(file, gz.ext)) | overwrite)) { 
		system(paste("pigz", file))
		return(paste0(basename(file), gz.ext))
	} else if (!grepl(paste0(gz.ext, "$"), file) & file.exists(paste0(file, gz.ext)) & rmIfExists) { 
		file.remove(file)
		return(paste0(basename(file), gz.ext))
	} else return(basename(file))
}

#check to see if file is gzipped and if not compress and return name
ugzipFile = function(file, overwrite = c(F), gz.ext = c(".gz")) {

	#Debug;
	#file = fqTemp1; overwrite = T; copy = F; gz.exts = ".gz";

	if (!is.null(file) && !is.na(file) && file != "") {

		outFile = gsub(gz.ext, "", file);

		if (file != outFile & (!file.exists(outFile) | overwrite)) { 
	#		system(paste0("gzip -", if (copy) "c", "d ", file, " > ", outFile))
			system(paste0("pigz -d ", file))
			return(outFile)
		} else print("Nothing")
	}
}
	
	
#check to see if fastq file is gzipped and if not compress and return name
gzipFastq = function(fqFile,  overwrite = c(F), rmIfExists = c(F)) {

	#Debug	fqFile = paste0(files$dir[i], files$fqMate1[i]);
	gz.ext = ".gz";
	fastq.ext = ".fastq"

	if (grepl(paste0(fastq.ext, "$"), fqFile) & (!file.exists(paste0(fqFile, gz.ext)) | overwrite)) { 
		system(paste("pigz", fqFile))
		return(paste0(basename(fqFile), gz.ext))
	} else if (grepl(paste0(fastq.ext, "$"), fqFile) & file.exists(paste0(fqFile, gz.ext)) & rmIfExists) { 
		file.remove(fqFile)
		return(paste0(basename(fqFile), gz.ext))
	} else return(basename(fqFile))
}

#QC file
fastqc = function(fqFiles, dir = c(NA), split = c(NA), threads = c(NA), outDir = c(NA)) {

	#Debug	fqFiles = files$fqMate1[i]; dir = files$dir[i]; split = "\\|"
	#fqFiles = paste0(files$dir[i], fqTemp1); dir = NA; split = NA; threads = 12;

	#require fqFiles
	if (!is.blank(fqFiles)) {

		#configure multiple fqFiles
		if (!is.na(split)) fqFiles = unlist(strsplit(fqFiles, split))
		if (!is.na(dir)) fqFiles = paste0(dir, fqFiles)
		if (length(fqFiles) > 1) fqFiles = paste(fqFiles, collapse = " ")	
	
		cmd = paste("fastqc", fqFiles)
		if (!is.na(threads)) cmd = paste(cmd , "-t", threads)
		if (!is.na(outDir)) {
			if (!file.exists(outDir)) dir.create(outDir)
			cmd = paste(cmd, "-o", outDir)
		}
		print(paste(cmd, Sys.time()))
		system(cmd)
	}
}

#function splits fq files to increase plexity of trimGalore call
splitFqs = function(fqFiles = c(NA), fqMates1 = c(NA), fqMates2 = c(NA), dir = c(NA), split = c("\\|"), lines = c(40000000)) {

	#Debug
	#fqFiles = files$fqFile[i]; fqMates1 = files$fqMate1[i]; fqMates2 = files$fqMate2[i]; dir = files$dir[i]; split = "\\|"; lines = 40000000;

	#set options so that lines is not put in scientific notation
	options("scipen" = 100)

	splitFqFiles = NULL
	splitFqMates1 = NULL
	splitFqMates2 = NULL

	splitCmd = paste("split -a 4 -d -l", lines)

	#If single-end reads
	if (!is.blank(fqFiles)) {

		#configure multiple fqFiles
		if (!is.na(split)) fqFiles = unlist(strsplit(fqFiles, split))

		for (fqFile in fqFiles) {
			print(paste("Splitting", fqFile, " - ", Sys.time()))

			#determine if gzipped
			if (grepl(".gz$", fqFile)) gz = T else gz = F

			if (gz) {
				system(paste0("zcat ", dir, fqFile, " | ", splitCmd, " - '", dir, gsub(".gz$", "-'", fqFile)))
				splitFqFiles = c(splitFqFiles, list.files(path = if (is.na(dir)) getwd() else dir, pattern = gsub(".gz$", "-", fqFile)))

			} else {
				system(paste0(splitCmd, " ", dir, fqFile, " '", dir, fqFile, "-'"))
				splitFqFiles = c(splitFqFiles, list.files(path = if (is.na(dir)) getwd() else dir, pattern = paste0(fqFile, "-")))
			}
		}

	} else if (!is.blank(fqMates1) & !is.blank(fqMates2)) {
		#configure multiple fqFiles
		if (!is.na(split)) {
			fqMates1 = unlist(strsplit(fqMates1, split))
			fqMates2 = unlist(strsplit(fqMates2, split))
		}

		for (fqMate1 in fqMates1) {
			print(paste("Splitting", fqMate1, " - ", Sys.time()))

			#determine if gzipped
			if (grepl(".gz$", fqMate1)) gz = T else gz = F

			if (gz) {
				system(paste0("zcat ", dir, fqMate1, " | ", splitCmd, " - '", dir, gsub(".gz$", "-'", fqMate1)))
				splitFqMates1 = c(splitFqMates1, list.files(path = if (is.na(dir)) getwd() else dir, pattern = gsub(".gz$", "-", fqMate1)))

			} else {
				system(paste0(splitCmd, " ", dir, fqMate1, " '", dir, fqMate1, "-'"))
				splitFqMates1 = c(splitFqMates1, list.files(path = if (is.na(dir)) getwd() else dir, pattern = paste0(fqMate1, "-")))
			}
		}

		for (fqMate2 in fqMates2) {
			print(paste("Splitting", fqMate2, " - ", Sys.time()))

			#determine if gzipped
			if (grepl(".gz$", fqMate2)) gz = T else gz = F

			if (gz) {
				system(paste0("zcat ", dir, fqMate2, " | ", splitCmd, " - '", dir, gsub(".gz$", "-'", fqMate2)))
				splitFqMates2 = c(splitFqMates2, list.files(path = if (is.na(dir)) getwd() else dir, pattern = gsub(".gz$", "-", fqMate2)))

			} else {
				system(paste0(splitCmd, " ", dir, fqMate2, " '", dir, fqMate2, "-'"))
				splitFqMates2 = c(splitFqMates2, list.files(path = if (is.na(dir)) getwd() else dir, pattern = paste0(fqMate2, "-")))
			}
		}
	}

	combine = gsub("^\\\\", "", split)
	return(list(paste0(splitFqFiles, collapse = combine), paste0(splitFqMates1, collapse = combine), paste0(splitFqMates2, collapse = combine)))
}

trimGalore = function(fqFiles = c(NA), fqMates1 = c(NA), fqMates2 = c(NA), dir = c(NA), split = c("\\|"), threads = c(NA), outDir = c(NA), qc = c(T), retainUnpaired = c(T), hardtrim5 = c(NA), hardtrim3 = c(NA), ...) { # 

	#Debug	fqFiles = files$fqFile[i]; fqMates1 = NA; fqMates2 = NA; dir = files$dir[i]; split = "\\|"; outDir = paste0(files$dir[i], "fqTrim/"); threads = 10; retainUnpaired = T;
	#Debug	fqFiles = fqFiles; fqMates1 = fqMates1; fqMates2 = fqMates2; dir = files$dir[i]; split = fqSplit; threads = threads; outDir = trimDir; retainUnpaired = T; fastqc = F;

	#require parallel package 
 	if (!is.na(threads) & !require("parallel")) stop()

	#determine if gzipped
	if (!is.blank(fqFiles) && grepl(".gz$", fqFiles)) gz = T else if (!is.blank(fqMates1) && grepl(".gz$", fqMates1)) gz = T else gz = F

	#determine the extension to make trimmed file names.  If the file ends with .fastq or .fq +/- .gz then trim_galore trims that off the file name, thus changing the naming convention 
	if (!is.blank(fqFiles)) {
		fqExt = gsub("^.*(\\.fq|\\.fastq)(\\.gz)?$", "\\1\\2", fqFiles)
		if (fqExt == fqFiles) trimFq = F else trimFq = T
	} else if (!is.blank(fqMates1)) {
		fqExt = gsub("^.*(\\.fq|\\.fastq)(\\.gz)?$", "\\1\\2", fqMates1)
		if (fqExt == fqMates1) trimFq = F else trimFq = T
	}
	
	if (is.na(hardtrim5) & is.na(hardtrim3)) {
		trimExt = paste0("_trimmed.fq", if (gz) ".gz")
		trim1Ext = paste0("_val_1.fq", if (gz) ".gz")
		trim2Ext = paste0("_val_2.fq", if (gz) ".gz")
		unpaired1Ext = paste0("_unpaired_1.fq", if (gz) ".gz")
		unpaired2Ext = paste0("_unpaired_2.fq", if (gz) ".gz")
	} else if (!is.na(hardtrim5)) {
		trimExt = paste0(".", hardtrim5, "bp_5prime.fq", if (gz) ".gz")
		trim1Ext = paste0(".", hardtrim5, "bp_5prime.fq", if (gz) ".gz")
		trim2Ext = paste0(".", hardtrim5, "bp_5prime.fq", if (gz) ".gz")
	} else if (!is.na(hardtrim3)) {
		trimExt = paste0(".", hardtrim3, "bp_3prime.fq", if (gz) ".gz")
		trim1Ext = paste0(".", hardtrim3, "bp_3prime.fq", if (gz) ".gz")
		trim2Ext = paste0(".", hardtrim3, "bp_3prime.fq", if (gz) ".gz")
	}

	#create out directory
	if (!missing(outDir) & !file.exists(outDir)) dir.create(outDir)

	#require fqFiles
	if (!is.blank(fqFiles)) {

		#configure multiple fqFiles
		if (!is.na(split)) fqFiles = unlist(strsplit(fqFiles, split))
		if (!is.na(dir)) fqFiles = paste0(dir, fqFiles)

		#for (fqFile in fqFiles) trimGaloreCall(fqFile = fqFile, ...)
		mcmapply(trimGaloreCall, fqFile = fqFiles, outDir = outDir, qc = qc, retainUnpaired = retainUnpaired, hardtrim5 = hardtrim5, hardtrim3 = hardtrim3, ..., mc.cores = threads)

		if (trimFq) fqFiles = gsub(fqExt, trimExt, fqFiles) else fqFiles = paste0(fqFiles, trimExt)
		if (!missing(outDir)) fqFiles = paste0(outDir, unlist(lapply(fqFiles, basename)))

		return(fqFiles)

	} else if (!is.blank(fqMates1) && !is.blank(fqMates2)) {

		#configure multiple fqMates
		if (!is.na(split)) {
			fqMates1 = unlist(strsplit(fqMates1, split))
			fqMates2 = unlist(strsplit(fqMates2, split))
		}
		if (!is.na(dir)) {
			fqMates1 = paste0(dir, fqMates1)
			fqMates2 = paste0(dir, fqMates2)
		}

		if (length(fqMates1) != length(fqMates2)) stop("Error: Different number of paired end files")

		#process trimGalore in parallel
		mcmapply(trimGaloreCall, fqMate1 = fqMates1, fqMate2 = fqMates2, outDir = outDir, qc = qc, retainUnpaired = retainUnpaired, hardtrim5 = hardtrim5, hardtrim3 = hardtrim3, ..., mc.cores = threads)

		#if the input fastq files has an extension recognized by trim Galore (i.e. .fastq or .fq +/- .gz) then trim that extension off the name
		if (trimFq) {
			if (is.na(hardtrim5) & is.na(hardtrim3)) fqUnpaired = c(gsub(fqExt, unpaired1Ext, fqMates1), gsub(fqExt, unpaired2Ext, fqMates2))
			fqMates1 = gsub(fqExt, trim1Ext, fqMates1)
			fqMates2 = gsub(fqExt, trim2Ext, fqMates2)
		} else {
			if (is.na(hardtrim5) & is.na(hardtrim3)) fqUnpaired = c(paste0(fqMates1, unpaired1Ext), paste0(fqMates2, unpaired2Ext))
			fqMates1 = paste0(fqMates1, trim1Ext)
			fqMates2 = paste0(fqMates2, trim2Ext)
		}
		
		if (!missing(outDir)) {
			fqMates1 = paste0(outDir, unlist(lapply(fqMates1, basename)))
			fqMates2 = paste0(outDir, unlist(lapply(fqMates2, basename)))
			if (is.na(hardtrim5) & is.na(hardtrim3)) fqUnpaired = paste0(outDir, unlist(lapply(fqUnpaired, basename)))
		}
		


		if (is.na(hardtrim5) & is.na(hardtrim3) & retainUnpaired) return(list(fqMates1, fqMates2, fqUnpaired)) else return(list(fqMates1, fqMates2)) 
	}
}


trimGaloreCall = function(fqFile = c(NA), fqMate1 = c(NA), fqMate2 = c(NA), outDir = c(NA), qc = c(T), retainUnpaired = c(T), hardtrim5 = c(NA), hardtrim3 = c(NA), ...) { # 

	#Debug	fqFile = fqFile; fqMate1 = NA; fqMate2 = NA; outDir = paste0(files$dir[i], "fqTrim/"); retainUnpaired = T

	cmd = "trim_galore"

	if (!is.na(hardtrim5)) cmd = paste(cmd, "--hardtrim5", hardtrim5)
	if (!is.na(hardtrim3)) cmd = paste(cmd, "--hardtrim3", hardtrim3)

	if (!is.na(outDir)) cmd = paste(cmd, "-o", outDir)
	if (qc) cmd = paste(cmd, "--fastqc")

	if (!is.na(fqFile)) {
		cmd = paste(cmd, fqFile)
	} else if (!is.na(fqMate1) & !is.na(fqMate2)) {
		cmd = paste(cmd, "--paired")
		if (retainUnpaired) cmd = paste(cmd, "--retain_unpaired")
		cmd = paste(cmd, fqMate1, fqMate2)
	}

	print(paste(cmd, Sys.time()))
	system(cmd)

}

#Sort bam file
sortBam = function(bamFile, bamSortFile = c(NA), delBam = c(F), threads = c(1), sortCmd = "samtools sort", mem = c("2G"), index = c(T)) {

	print(paste("sorting Bam", basename(bamFile), Sys.time()))

	idxExt = ".bai";
	bamExt = ".bam";
	bamSortExt = ".sort.bam";

	if (is.na(bamSortFile)) bamSortFile = gsub(bamExt, bamSortExt, bamFile)

	#sort BAM file
	sortCmd = paste(sortCmd, "-m", mem, "-@", threads, "-o", bamSortFile, bamFile)
	system(sortCmd)
	
	#Remove unsorted bam file
	if (delBam) file.remove(bamFile)
	if (delBam & file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))

	#Create BAM index for fast access
	if (index) indexBam(bamSortFile)		

	return(bamSortFile)
	
}

#Use samtools to index bam file 
indexBam = function(bamFile, threads = c(1)) {

	idxCmd = paste("samtools index -@", threads, bamFile)
	system(idxCmd)
}

#Use samtools to mark duplicates in bam file 
markDups = function(bamFile, delBam = c(F), threads = c(1)) {
	#Debug:
	#bamFile = paste0(files$dir[i], files$bamFile[i]); delBam = T; threads = threads;

	idxExt = ".bai";
	bamExt = ".bam";
	bamDupExt = ".dupMark.bam";
	bamFixExt = ".fix.bam";
	
	bamFixFile = gsub(paste0(bamExt, "$"), bamFixExt, bamFile)
	bamDupFile = gsub(paste0(bamExt, "$"), bamDupExt, bamFile)
	
	bamSortFile = sortBam(bamFile, threads = threads, sortCmd = "samtools sort -n", index = F)
	
	print(paste("Fixing", bamFile, Sys.time()))

	fixCmd = paste("samtools fixmate -@", threads, "-m", bamSortFile, bamFixFile)
	system(fixCmd)
	
	bamSortFixFile = sortBam(bamFixFile, threads = threads)

	print(paste("Marking Duplicates", bamFile, Sys.time()))
	dupCmd = paste0("samtools markdup -s -@ ", threads, " ", bamSortFixFile, " ", bamDupFile)
	system(dupCmd)

	#Create BAM index for fast access
	indexBam(bamDupFile)		

	#Remove BAM without duplicates marked
	if (delBam) {
		file.remove(bamFile)
		file.remove(bamSortFile)
		file.remove(bamFixFile)
		file.remove(bamSortFixFile)

		if (file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))
		if (file.exists(paste0(bamSortFile, idxExt))) file.remove(paste0(bamSortFile, idxExt))
		if (file.exists(paste0(bamFixFile, idxExt))) file.remove(paste0(bamFixFile, idxExt))
		if (file.exists(paste0(bamSortFixFile, idxExt))) file.remove(paste0(bamSortFixFile, idxExt))
	}

	return(basename(bamDupFile))
}

#Use picard to mark duplicates in bam file 
markDupsPICARD = function(bamFile, picardCmd, delBam = c(F)) {

	#Debug:
	#bamFile = paste0(files$dir[i], files$bamFile[i]); picardCmd = picardCmd; delBam = T;

	print(paste("Marking Duplicates", bamFile, Sys.time()))

	idxExt = ".bai";
	bamExt = ".bam";
	bamDupMarkExt = ".dupMark.bam";
			
	bamDupFile = gsub(paste0(bamExt, "$"), bamDupMarkExt, bamFile)

	#Mark duplicates using picard
	system(paste0(picardCmd, "INPUT=", bamFile, " OUTPUT=", bamDupFile, " METRICS_FILE=", bamDupFile, ".metrics"))

	#Create BAM index for fast access
	indexBam(bamDupFile)		

	#Remove BAM without duplicates marked
	if (delBam) file.remove(bamFile)
	if (delBam & file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))

	return(bamDupFile)
}


#Get read counts for bam file
getBamCts = function(bamFile, baiExt = c(".bai")) {
	#Debug;
	#bamFile = paste0(files$dir[i], files$bamFile[i]); unMappedBamFile = NA; baiExt = ".bai"; 

	print(paste("Getting Bam File Read Counts", basename(bamFile), Sys.time()))

	#create index for bamFile if it doesn't exist
	if (!file.exists(paste0(bamFile, baiExt))) indexBam(bamFile); 

	#total reads (first + mate)
	stat = idxstatsBam(bamFile)

	return(c(sum(stat$mapped), sum(stat$unmapped)))
}

#Function that reads Bismark Call table and returns data
readBismarkCall = function(file, gz = c(NA)) {

	#read in ot and ob strand sort strand CpG calls
	covNames = c("read", "fauxStrand", "chr", "start", "meth")
	covClasses = c("NULL", "NULL", "character", "numeric", "character")

	if ((is.na(gz) & grepl("\\.gz$", file)) | gz) call = gzfile(file, "r") else call = file(file, "r")
	data = suppressWarnings(read.table(call, skip = 1, header = F, sep = "\t", comment.char = "", colClasses = covClasses, col.names = covNames))
	close(call)

	return(data)
}

#Function that reads Bistools table and returns data
readBistoolsCall = function(file, gz = c(NA), covClasses = c(NA), threads = c(10)) {

	#file = "/home/bbarwick/Documents/Boss/Dnmt3ab/RRBS/RRBS_fasta/Dnmt3_702/t.call.gz"

	#read in ot and ob strand sort strand CpG calls
	covNames = c("readPos", "flag", "rname", "strand", "pos", "qwidth", "cigar", "meth")
	if (all(is.na(covClasses))) covClasses = c("NULL", "NULL", "character", "NULL", "numeric", "NULL", "NULL", "integer")

	if ((is.na(gz) & grepl("\\.gz$", file)) | gz) data = fread(cmd = paste("zcat <", file), colClasses = covClasses, nThread = threads) else data = fread(file, colClasses = covClasses, nThread = threads)
	gc()
	return(data)
}

#Function that reads Bistools Coverage data table
readBistoolsCov = function(file, gz = c(NA), covClasses = c(NA)) {

	#Debug
	#file = callFile; covClasses = NA; gz = NA;

	#read in ot and ob strand sort strand CpG calls
	if (any(is.na(covClasses))) covClasses = c("character", "numeric", "numeric", "integer", "integer")

	if ((is.na(gz) & grepl("\\.gz$", file)) | gz) cov = fread(cmd = paste("zcat <", file), colClasses = covClasses) else cov = fread(file, colClasses = covClasses)
	setkey(cov, rname, pos)

	return(cov)
}

#Function that writes Bistools Coverage to tab delimited txt
writeBistoolsCov = function(cov, file, compress = c(NA), overwrite = c(T), sep = c("\t"), round = options("digits")[[1]], ...) {

	#Debug
	#cov = sum; file = paste0(covDir, grp, ".cov.", minCov, ".cov"); gz = NA; sep = "\t"

	if (is.blank(compress) & grepl("\\.gz$", file)) {
		compress = T
		file = gsub("\\.gz$", "", file)
	} else compress = F

	if (!is.na(round) & is.integer(round)) {
		numericCols = sapply(cov, class)
		numericCols = names(numericCols[numericCols == "numeric"])
		cov[, (numericCols):=round(.SD, round), .SDcols = numericCols]
	}

	fwrite(cov, file = file, sep = sep, row.names = F, ...)
	if (compress) system(paste("pigz", if (overwrite) "-f", file))
}

#Function that writes Bistools Coverage to tab delimited bed file
writeBistoolsCovToBed = function(cov, file, compress = c(T), sep = c("\t"), overwrite = c(T), chrCol = c("rname"), posCol = c("pos")) {

	#Debug
	#cov = sum; file = paste0(covDir, grp, ".cov.", minCov, ".cov"); gz = NA; sep = "\t"	
	
	#handle gz compression
	if ((is.na(compress) & grepl("\\.gz$", file)) | (!is.na(compress) & compress)) {
		compress = T
		file = gsub("\\.gz$", "", file)
	}

	#options(scipen = 999); #no scientific notation
	covBed = data.table(chr = cov[[chrCol]], start = format(cov[[posCol]] - 1, scientific = F), end = format(cov[[posCol]], scientific = F)); #, strand = rep("*", dim(cov)[1])

	fwrite(covBed, file = file, sep = sep, col.names = F, row.names = F, quote = F)

	#handle gz compression
	if (compress) system(paste("pigz", if (overwrite) "-f", file))

}

writeBistoolsCovToBw = function(cov, file, chrCol = c("chr"), posCol = c("pos"), scoreCol = c("meth"), removeChrs = c(NA), removeBg = c(T), digits = c(options("digits")[[1]]), ...) {

	#Debug
	#cov = cov; file = paste0(covDir, files$sample[i], ".minCov", minCov, ".minSamples", minSampleCov, ".bw"); chrCol = "rname"; posCol = "pos"; scoreCol = paste0(files$sample[i], ".meth"); genome = genome; removeChrs = removeChrs; removeBg = T; digits = 3
	
	bgExt = ".bedGraph";

	covBg = data.table(chr = cov[[chrCol]], start = cov[[posCol]] - 1, end = cov[[posCol]], score = cov[[scoreCol]])
	if (!is.na(removeChrs)) covBg = subset(covBg, !grepl(removeChrs, covBg$chr))
	covBg = subset(covBg, !is.na(covBg$score))
	covBg = covBg[, lapply(.SD, format, scientific = F)]
	covBg$score = round(as.numeric(covBg$score), digits)

	fwrite(covBg, file = paste0(file, bgExt), sep = "\t", col.names = F, row.names = F, quote = F)

	bedgraph2bw(paste0(file, bgExt), bwFile = file, ...);
	
	if (removeBg) file.remove(paste0(file, bgExt))
}

#Converts a BisTools Cov data.table to bedGraph
bistoolsCov2Bedgraph = function(bg, bgFile = c(NA), minCov = c(1), removeChrs = c(NA)) {

	#bg = cov; bgFile = "test.bedGraph.gz"; minCov = 10; removeChrs = c("gi|215104|gb|J02459.1|LAMCG")

	#Calculate total coverage and filter for minCov
	bg[, start:=pos-1]
	bg[, end:=pos]
	bg = bg[m+u >= minCov]
	
	#remove extraneous columns
	bg[, m:=NULL]
	bg[, u:=NULL] 
	bg[, pos:=NULL]

	#reorder columns
	setcolorder(bg, c("rname", "start", "end", "meth"))

	#remove chromsomes contained in the removeChrs variable
	if (any(!is.na(removeChrs))) bg = bg[!(rname %in% removeChrs)]

	#sort by chromosome and position
	setkeyv(bg, c("rname", "start"))
	
	#if bgFile is given write to bgFile else return bg
	if (!is.na(bgFile)) {
		options(scipen=999)
		if (grepl("\\.gz$", bgFile)) bgCon = gzfile(bgFile, open = "w") else bgCon = file(bgFile, open = "w")
		write.table(bg, file = bgCon, sep = "\t", row.names = F, quote = F, col.names = F)
		close(bgCon)
	} else return(bg)
}

##Wrapper function to convert bedgraph to bigwig
bedgraph2bw = function(bgFile, genome = c(NA), si = c(NA), chromFile = c(NA), bwFile = c(gsub(".bedGraph.*$", ".bw", bgFile))) {

	#bgFile = paste0(files$dir[i], files$bgFile[i]); genome = genome; seqLengths = NA; chromFile = NA; bwFile = gsub(".bedGraph.*$", ".bw", bgFile)	
	#bgFile = paste0(file, bgExt); bwFile = file; genome = genome; si = NA; chromFile = NA
	
	#write chromosome sizes to txt file for bedGraphToBigWig 
	if (suppressWarnings(!is.na(genome))) {
		chromFile = paste0(genome@pkgname, ".chrom.txt")
		write.table(seqlengths(genome), file = chromFile, sep = "\t", quote = F, col.names = F)
	} else if (!is.na(suppressWarnings(si))) {
		if (is.na(chromFile)) chromFile = paste0(si@genome[1], ".chrom.txt")
		write.table(cbind(si@seqnames, si@seqlengths), file = chromFile, sep = "\t", quote = F, col.names = F, row.names = F);
	}

	if (!is.na(chromFile)) {
		bwCmd = paste("bedGraphToBigWig", gsub(".gz$", "", bgFile), chromFile, bwFile)
		system(bwCmd)
	} else warning("No compatible chromosome sequence length data");
	
}

#Wrapper function to convert bedgraph to bigBed
bed2bb = function(bedFile, genome, bbFile = c(gsub(".bed.*$", ".bb", bedFile))) {

	#make bigwig file. decompress bedgraph if necessary
	if (grepl(".gz$", bedFile)) system(paste("gzip -cd", bedFile, ">", gsub(".gz$", "", bedFile)))	

	#write chromosome sizes to txt file for bedGraphToBigWig 
	write.table(seqlengths(genome), file = paste0(genome@pkgname, ".chrom.txt"), sep = "\t", quote = F, col.names = F)

	bbCmd = paste("bedToBigBed", gsub(".gz$", "", bedFile), paste0(genome@pkgname, ".chrom.txt"), bbFile)
	system(bbCmd)
	
	if (grepl(".gz$", bedFile)) file.remove(gsub(".gz$", "", bedFile))
}

##### Wrapper function for bismarkBam2call specifically for CpGs
bismarkBam2callCpg = function(...) bismarkBam2call(callChar = c("z"), ...)
bismarkBam2callChh = function(...) bismarkBam2call(callChar = c("x"), ...)
bismarkBam2callChg = function(...) bismarkBam2call(callChar = c("h"), ...)

##### Helper function which checks to see if bamFile is indexed and if so reads in bam file 1 chromosome at a time converting it to a methylation calls data table
bismarkBam2call = function(bamFile, callFile = c(NA), covFile = c(NA), collapseCpgFile = c(NA), collapseCpgCovFile = c(NA), trim1 = c(0, 0), trim2 = c(0, 0), trimFile = c(NA), callChar = c("z"), compress = c(T), threads = c(10), bWhat = c("rname", "strand", "pos", "flag", "qwidth", "cigar"), maxRegionSize = c(Inf), overwrite = c(F), ... ) {
	#Debug
	#bamFile = paste0(files$dir[i], files$bamFile[i]); callFile = NA; covFile = NA; collapseCpgFile = paste0(files$dir[i], files$collapseCpgFile[i]); collapseCpgCovFile = paste0(files$dir[i], files$covFile[i]); trim1 = c(0, 0); trim2 = c(0, 0); trimFile = NA; callChar = "z"; compress = T; threads = threads; bWhat = c("rname", "strand", "pos", "flag", "qwidth", "cigar"); maxRegionSize = 1e7

	bf = BamFile(bamFile)
	bParam = ScanBamParam(what = bWhat, tag = "XM")

	if (length(index(bf)) == 0 | is.na(index(bf))) {
		print(paste(bamFile, "is not indexed reading in entire bam file ..."));
		call = bismarkBam2callDt(bamFile, callFile = callFile, collapseCpgFile = collapseCpgFile, trim1 = trim1, trim2 = trim2, trimFile = trimFile, callChar = callChar, bParam = bParam, adj = adj, col.names = T); #compress = compress, threads = threads,
	} else {
		print(paste(bamFile, "is indexed reading in one chromosome at a time ..."));
		chrLen = scanBamHeader(bf)$targets
		chrs = names(chrLen)
		
		if (!is.na(callFile)) chrCallFiles = paste0(make.names(chrs), ".", basename(callFile)) else chrCallFiles = rep(NA, length(chrs))
		if (!is.na(covFile)) chrCovFiles = paste0(make.names(chrs), ".", basename(covFile)) else chrCovFiles = rep(NA, length(chrs))
		if (!is.na(collapseCpgFile)) chrCollapseCpgFiles = paste0(make.names(chrs), ".", basename(collapseCpgFile)) else chrCollapseCpgFiles = rep(NA, length(chrs))
		if (!is.na(collapseCpgCovFile)) chrCollapseCpgCovFiles = paste0(make.names(chrs), ".", basename(collapseCpgCovFile)) else chrCollapseCpgCovFiles = rep(NA, length(chrs))
		if (!is.na(trimFile)) chrTrimFile = paste0(make.names(chrs), ".", basename(trimFile)) else chrTrimFiles = rep(NA, length(chrs))

		####
		col.names = c(T, rep(F, length(chrs) - 1))
		bParams = list()
		for (i in 1:length(chrs)) {
			bamWhich(bParam) = GRanges(seqnames = chrs[i], ranges = IRanges(start = 0, end = chrLen[which(names(chrLen) == chrs[i])]))			
			bParams[[i]] = bParam
		}

		#parallel call to bismarkBam2callDt
		#mcmapply(bismarkBam2callDt, bamFile = rep(bamFile, length(chrs)), callFile. = chrCallFiles, collapseCpgFile. = chrCollapseCpgFiles, collapseCpgCovFile. = chrCollapseCpgCovFiles, trim1 = rep(trim1, length(chrs)), trim2 = rep(trim2, length(chrs)), trimFile = chrTrimFiles, callChar. = callChar, bParam = bParams, maxRegionSize = rep(maxRegionSize, length(chrs)), col.name = col.names, message = chrs, ..., mc.cores = threads)
		mcmapply(bismarkBam2callDt, bamFile = bamFile, callFile. = chrCallFiles, covFile. = chrCovFiles, collapseCpgFile. = chrCollapseCpgFiles, collapseCpgCovFile. = chrCollapseCpgCovFiles, trim1 = trim1, trim2 = trim2, trimFile = chrTrimFiles, callChar. = callChar, bParam = bParams, maxRegionSize = maxRegionSize, col.name = col.names, message = chrs, ..., mc.cores = threads)

		#combine chr files
		if (!is.na(callFile)) catFiles(paste0(chrCallFiles[file.exists(chrCallFiles)], collapse = "|"), outFile = callFile, compress = T, del = T, overwrite = overwrite)
		if (!is.na(collapseCpgFile)) catFiles(paste0(chrCollapseCpgFiles[file.exists(chrCollapseCpgFiles)], collapse = "|"), outFile = collapseCpgFile, compress = T, del = T, overwrite = overwrite)
		if (!is.na(collapseCpgCovFile)) catFiles(paste0(chrCollapseCpgCovFiles[file.exists(chrCollapseCpgCovFiles)], collapse = "|"), outFile = collapseCpgCovFile, compress = T, del = T, overwrite = overwrite)
		if (!is.na(trimFile)) catFiles(paste0(chrTrimFiles[file.exists(chrTrimFiles)], collapse = "|"), outFile = trimFile, del = T, overwrite = overwrite)
	}
}

##### Helper function for converting a bamFile to methylation calls data table
#Alternative default bParam = c(ScanBamParam(what = c("qname", "rname", "strand", "pos", "qwidth"), tag = "XM"))

bismarkBam2callDt = function(bamFile, callFile. = c(NA), covFile. = c(NA), collapseCpgFile. = c(NA), collapseCpgCovFile. = c(NA), trim1 = c(0, 0), trim2 = c(0, 0), trimFile = c(NA), callChar. = c("z"), bParam = c(ScanBamParam(what = c("rname", "strand", "pos", "qwidth", "cigar"), tag = "XM")), maxRegionSize = c(Inf), scanBamBuffer = c(1e4), col.names = c(F), message, logFile = c(NA), ... ) {
	#bamFile = bamFile; callFile. = chrCallFiles[i]; covFile. = NA; collapseCpgFile. = chrCollapseCpgFiles[i]; collapseCpgCovFile. = chrCollapseCpgCovFiles[i]; trim1 = trim1; trim2 = trim2; trimFile = chrTrimFiles[i]; callChar. = "z"; bParam = bParams[[i]]; maxRegionSize = maxRegionSize; scanBamBuffer = 1e4; col.names = col.names[i]; message = chrs[i];
	
	if (!is.na(message)) cat(paste(message, "...", Sys.time(), "\n"))

	#set chrom parameters
	start = start(bamWhich(bParam)[[1]])
	end = end(bamWhich(bParam)[[1]])
	
	if (width(bamWhich(bParam)[[1]]) > maxRegionSize) {
		starts = as.integer(seq(from = start, to = end, by = maxRegionSize))
		ends = as.integer(c(seq(from = start + maxRegionSize - 1, to = end, by = maxRegionSize), end))
	} else {
		starts = start
		ends = end
	}		

	for (i in 1:length(starts)) {
	
		start(bamWhich(bParam)[[1]]) = max(c(starts[i] - scanBamBuffer, 0)) # add scanBamBuffer range to account for insertions / deletions overlapping range
		end(bamWhich(bParam)[[1]]) = min(c(ends[i] + scanBamBuffer, end))

		if (!is.na(logFile)) memLog(paste("bismarkBam2callDt scanning Bam region", message, starts[i], "-", ends[i]), file = logFile)

		#Read in from bam file
		bam = scanBam(bamFile, param = bParam)[[1]]

		### Get CpG methylation Calls
		m = strsplit(gsub(paste0("[^", callChar., "]"), "", bam$tag$XM, ignore.case = T), ""); #Parse CpG Methylation Calls "Z" and "z" from bam file
	
		#Ensure there are records prior and methylation calls
		if (length(bam[[1]]) > 0 & !all(lapply(m, length) == 0)) {
		
			mLoc = gregexpr(callChar., sequenceLayer(BStringSet(bam$tag$XM), bam$cigar), ignore.case = T); #Parse C methylation calls and save locations of Z and z to list of lists. Use sequenceLayer to account for insertions and deletions
			mLoc = mLoc[unlist(lapply(mLoc, function(x) return(!all(x == -1))))]; #Remove records with no call (i.e. neither z or Z) for which gregexpr returns -1
			call = data.table(mCall = unlist(m), readPos = unlist(mLoc))

			#add all other field in bamWhat 
			fields = names(bam)[names(bam) != "tag"]
			for (field in fields) call[, eval(parse(text = paste0(field, ":=rep(bam[[field]], lapply(m, length))")))]
			rm(m, mLoc); gc();

			#Identify second mate pair reads for position adjustment
			call[, secondMate := bamFlagTest(flag, "isSecondMateRead")]

			#adjust position by read position 
			call[, pos:=pos+readPos-1]

			#adjust second mate reads by 1 strand-specific
			call[strand=="+" & secondMate, pos:=pos-1]; 
			call[strand=="-" & secondMate, pos:=pos+1];

			#order by pos
			call = call[order(pos)]

			#Convert mCall (e.g. "Z" and "z") into boolean (i.e. 1 and 0) 
			call[, meth:=0]
			call[mCall==toupper(callChar.), meth:=1]
			call[, mCall:=NULL];		#eliminate unnecessary columns
			gc()

			if (!is.na(logFile)) memLog(paste("bismarkBam2callDt extracted meth calls region", i), file = logFile)

			#trim left and right sides of read, mate specific
			if (length(trim1) == 2 && !all(trim1 == 0)) {
				mtrim1 = getCallDtTrimmings(subset(call, !secondMate), trim = trim1)
				call = rbind(trimCallDt(subset(call, !secondMate), trim = trim1), subset(call, secondMate))
				print("trim1")
			} else mtrim1 = NULL

			if (length(trim2) == 2 && !all(trim2 == 0)) {
				mtrim2 = getCallDtTrimmings(subset(call, secondMate), trim = trim2)
				call = rbind(subset(call, !secondMate), trimCallDt(subset(call, secondMate), trim = trim2))
				print("trim2")
			} else mtrim2 = NULL
			
			if (!is.null(mtrim1) & !is.null(mtrim2)) { mtrim = rbind(mtrim1, mtrim2); print("trim") } else mtrim = NULL

			if (i == 1) append = F else { col.names = F; append = T }

			if (!is.na(trimFile) & !is.null(mtrim)) fwrite(mtrim, file = trimFile, sep = "\t", row.names = F, quote = F, col.names = col.names, append = append)
			if (!is.na(callFile.)) {
				fwrite(subset(call, pos >= starts[i] & pos <= ends[i]), file = callFile., sep = "\t", row.names = F, quote = F, col.names = col.names, append = append); #write out call file, remove positions not in start/end range 
				if (!is.na(covFile.)) bistoolsCall2Cov(call, covFile., col.names = col.names, append = append, logFile = logFile)
			}
	
			#write out call file, remove positions not in start/end range
			if (!is.na(collapseCpgFile.)) {
				call$pos[call$strand == "-"] = call$pos[call$strand == "-"] - 1
				call = subset(call, pos >= starts[i] & pos <= ends[i]); #

				fwrite(call, file = collapseCpgFile., sep = "\t", row.names = F, quote = F, col.names = col.names, append = append)

				if (!is.na(collapseCpgCovFile.)) bistoolsCall2Cov(call, collapseCpgCovFile., col.names = col.names, append = append, logFile = logFile)
			}

			suppressWarnings(rm(call, mtrim1, mtrim2, mtrim)); gc()
		}			
	}
}

#Summarizes Bismark Call data to coverage
bistoolsCall2Cov = function(call, covFile = c(NA), col.names = c(T), append = c(F), logFile = c(NA)) {

	#Debug
	#call = readBistoolsCall(paste0(files$dir[i], files$collapseCpgFile[i])); covFile = paste0(files$dir[i], files$covFile[i]);

	if (!is.na(logFile)) memLog("bistoolsCall2Cov start", file = logFile)	

	#use data.table to summarize average methylation 'meth', and total calls of methylated 'm' and 'u' states	
	mean = data.table(call)[, list(meth = mean(as.numeric(meth))), by = "rname,pos"]
	u = data.table(call)[meth == 0, list(u = length(meth)), by = "rname,pos"]
	m = data.table(call)[meth == 1, list(m = length(meth)), by = "rname,pos"]

	rm(call); gc(); #free memory

	if (!is.na(logFile)) memLog("bistoolsCall2Cov mean, u, m calculated", file = logFile)	
	
	#merge data
	cov = merge(mean, u, by = c("rname", "pos"), all = T)
	cov = merge(cov, m, by = c("rname", "pos"), all = T)

	rm(mean); rm(u); rm(m); gc(); #free memory

	#replace NA values with 0 for the m and u counts
	cov$u[is.na(cov$u)] = 0
	cov$m[is.na(cov$m)] = 0

	if (!is.na(logFile)) memLog("bistoolsCall2Cov cov table made", file = logFile)	

	#if bgFile is given write to bgFile else return bg
	if (!is.na(covFile)) writeBistoolsCov(cov, covFile, col.names = col.names, append = append) else invisible(cov)
	rm(cov); gc(); #free memory

}

#Function to check if bam file is paired
bamPaired = function(bamFile, yieldSize = c(1000)) {
	#Debug;
	#bamFile = bamFile; yieldSize = 1000;

	bam = scanBam(BamFile(bamFile, yieldSize = yieldSize))[[1]]
	sum = colSums(bamFlagAsBitMatrix(bam$flag))
	if (sum[["isPaired"]] > 0) return(T) else return(F)
}

##### Helper function for directionally trimming a methylation calls data table 
trimCallDt = function(call, trim = c(0, 0)) {

	#Debug
	#call = m1call; trim = trim1;

	call = call[((readPos > trim[1] & strand == "+") | (qwidth - readPos >= trim[1] & strand == "-")) & ((qwidth - readPos >= trim[2] & strand == "+") | (readPos > trim[2] & strand == "-"))]
	return(call)
}

##### Helper function for getting directionally trimmed a methylation calls data table 
getCallDtTrimmings = function(call, trim = c(0, 0)) {

	#Debug
	#call = m1call; trim = trim1;
	t1 = call[(readPos <= trim[1] & strand == "+") | (qwidth - readPos < trim[1] & strand == "-")]
	t2 = call[(qwidth - readPos < trim[2] & strand == "+") | (readPos <= trim[2] & strand == "-")]

	return(rbind(t1, t2))
}

#Function gets estimated fragment sizes by chromosome for a bam file and weights it by the # of reads / chromosome
getBamFragSizes = function(bamFiles, sbf = c(scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F)), fragMethod = c("SISSR"), threads = c(1), removeChrs = c(NA)) {
	#Debug	
	#bamFiles = paste0(files$dir[i], files$bamFile[i]); sbf = scanBamFlag(isSecondaryAlignment = F, isDuplicate = F); fragMethod = "SISSR"; threads = 6;
	
	fragSizes = NULL;

	for (i in 1:length(bamFiles)) {

		chrs = unique(getBamChrs(bamFiles[i]))		
		if (!is.na(removeChrs)) chrs = chrs[!grepl(removeChrs, chrs)]
		si = seqinfo(BamFile(bamFiles[i]));

		#set up variables to pass parallel call of getBamFragSize
		messages = paste(basename(bamFiles[i]), chrs)
		sbps = list()
		for (chr in chrs) sbps[[chr]] = ScanBamParam(flag = sbf, which = GRanges(seqnames = chr, ranges = IRanges(start = 1, end = si@seqlengths[si@seqnames == chr])))

		#estimate fragment size if bam is not paired
		if (!bamPaired(bamFiles[i])) {
			fragSizeChr = mcmapply(getBamFragSize, bamFile = bamFiles[i], sbp = sbps, fragMethod = fragMethod, message = messages, USE.NAMES = F, mc.cores = threads); #parallel call to get getBamFragSizes
		} else {
			fragSizeChr = mcmapply(getBamISize, bamFile = bamFiles[i], sbp = sbps, message = messages, USE.NAMES = F, mc.cores = threads); #parallel call to get getBamISize
		}

		#calc average fragsize weighted by chromosomal reads and add that to the list fragSizes
		fragSize = sum(as.numeric(names(fragSizeChr)) * fragSizeChr, na.rm = T) / sum(as.numeric(names(fragSizeChr)))
		fragSizes = c(fragSizes, fragSize);
	}
	return(fragSizes)
}

#function to get Bam insert size
getBamISize = function(bamFile, sbp = c(ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F))), maxISize = c(2000), message = c(bamFile), verbose = c(T)) {

	#Debug
	#bamFile = bamFiles[i]; sbp = sbps[[1]]; maxISize = 2000; message = messages[[1]];

	#import insert size
	sbp@what = c(sbp@what, "isize")
	bam = scanBam(bamFile, param = sbp)[[1]]

	if (length(bam) > 0) fragSize = mean(abs(bam$isize)[abs(bam$isize) <= maxISize], na.rm = T) else fragSize = NA	
	names(fragSize) = length(bam$isize)

	if (verbose) print(paste(message, "length: ", length(bam), "insert size: ", fragSize, Sys.time()));

	return(fragSize)
}



#function to get Bam fragment size used 
getBamFragSize = function(bamFile, sbp = c(ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F))), fragMethod = c("SISSR"), message = c(bamFile), verbose = c(T)) {
	#Debug
	#bamFile = bamFiles[i]; sbp = sbps[[1]]; fragMethod = "SISSR"; message = messages[[1]];

	bam = granges(readGAlignments(bamFile, param = sbp))

	if (length(bam) > 0) fragSize = estimate.mean.fraglen(bam, method = fragMethod) else fragSize = NA
	names(fragSize) = as.character(length(bam))

	if (verbose) print(paste(message, "length: ", length(bam), "estimated fragment size: ", fragSize, Sys.time()));

	return(fragSize)
}



#Function to return chromosomes included in a Bam file
getBamChrs = function(bamFile) return(names(scanBamHeader(bamFile)[[1]][1]$targets));

getBamChrLengths = function(bamFile) {

	#Debug: bamFile = "/home/bbarwick/Documents/Boss/Bcell/MeDIPseq/Bcell1/Bcell1.MeDIP.unique.bam";

	h = scanBamHeader(bamFile)
	cl = h[[1]][2]$text[names(h[[1]][2]$text) == "@SQ"]
	
	chrs = gsub("SN:", "", unlist(cl)[grepl("SN:", unlist(cl))]);
	len = as.numeric(gsub("LN:", "", unlist(cl)[grepl("LN:", unlist(cl))]));
	names(len) = chrs;

	return(len);
}


## bamToBigWig constants
gFragSize = NA;
gBwFile = "gsub('.bam', paste('.frag', paste(round(fragSizes, 0), collapse = '_'), '.norm', format(normReads, scientific = FALSE), '.bw', sep = ''), bamFiles)";

#bamToBigWig function: converts a bam file into a bigWig coverage file.  Removes chromosomes that match the variable removeChrs by grepl.  Normalizes for read count .Only reads not in removeChrs are counted. The normReads variable sets the sequencing coverage. The default normalization is reads per million (e.g. normRead = 1e6).  
bamToBigWig = function(bamFiles, bwFile = c(gBwFile), fragSizes = c(gFragSize), sigDigits = c(gSigDigits), readCounts = c(NA), normReads = c(gNormReads), extCall = c(T), sbf = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F), removeChrs = c(NA), threads = c(1), ...) {
	#Debug
	#bamFiles = paste0(files$dir[i], files$bamFile[i]); bwFile = gBwFile; fragSizes = files$fragSize[i]; sigDigits = gSigDigits; normReads = 1e6; normFactor = NA; extCall = T; sbf = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F); removeChrs = removeChrs; threads = threads
	#bamFiles = paste0(files$dir[i], files$bamFile[i]); bwFile = paste0(files$dir[i], files$sample[i], ".rpm.bw"); fragSizes = gFragSize; sigDigits = gSigDigits; readCounts = NA; normReads = gNormReads; normFactor = NA; flag = sbf; removeChrs = "gl00|chrUn"; extCall = T
	#bamFiles = bamFiles; bwFile = paste0(covDir, grp, ".union.rppm.bw"); fragSizes = gFragSize; sigDigits = gSigDigits; readCounts = sum(filesGrp$unique.incChr.reads); normReads = gNormReads / frip; sbf = sbf; threads = threads; removeChrs = removeChrs;	
	#bamFiles = paste0(files$dir[i], files$bamFile[i]); bwFile = paste0(outDir, make.names(files$sample[i]), ".union.rppm.bw"); fragSizes = gFragSize; sigDigits = gSigDigits; readCounts = files$nondup.chrInc.reads[i]; normReads = gNormReads / files[[paste0("frip.", sig)]][i]; sbf = flag; removeChrs = removeChr; threads = threads;

	#if the number of bamFiles & fragSizes aren't equal then through warning
	if (!all(is.na(fragSizes)) & length(bamFiles) != length(fragSizes)) warning("Number of bamFiles != number of fragSizes")
 
	#Get genome version and seqinfo
	si = list(); for (i in 1:length(bamFiles)) si[[i]] = seqinfo(BamFile(bamFiles[i]))
	chrs = unique(unlist(lapply(si, seqnames))); 	#get all chromosomes in bamFiles
	chrs = chrs[order(chrs)]
	if (!is.na(removeChrs)) chrs = chrs[!grepl(removeChrs, chrs)]; 	#remove chromosomes that match removeChrs

	gr = GRanges(seqnames = chrs,  ranges = IRanges(start = rep(0, length(chrs)), end = si[[1]]@seqlengths[match(chrs, si[[1]]@seqnames)]))

	#Get readcounts for each bamFile
	if (all(is.na(readCounts))) {
		readCounts = list()
		for (i in 1:length(bamFiles)) readCounts[[i]] = sum(countBam(bamFiles[i], param = ScanBamParam(flag = sbf, which = gr))$records)
	}

	#Create bigwig and bedfile output file names	
	if (bwFile == gBwFile) bwFile = eval(parse(text = gBwFile))
	bgFile = paste0(bwFile, ".bedGraph")
	bgFiles = paste0(paste0(bwFile, ".bedGraph"), ".", chrs)
	
	#set up scan bam params
	sbps = list()	
	for (chr in chrs) sbps[[chr]] = ScanBamParam(flag = sbf, which = gr[which(chrs == chr)])

	print(paste("bamToBedgraph", Sys.time()))

	#parallel call to bismarkBam2callDt
	mcmapply(bamToBedgraph, bamFiles. = rep(list(bamFiles), length(chrs)), bgFile. = bgFiles, fragSizes. = fragSizes, sigDigits. = sigDigits, normFactors. = rep(list(normReads / unlist(readCounts)), length(chrs)), sbp. = sbps, message. = chrs, mc.cores = threads) 

	#concatenate and clean-up
	catFiles(bgFiles, bgFile)
	file.remove(bgFiles)

	#if external calll write to bedgraph then call Kent tools
	if(extCall) {
		tmpFile = paste0(bwFile, ".temp")
		writeSiChromSizes(si[[1]], tmpFile, removeChrs)
		system(paste("bedGraphToBigWig", bgFile, tmpFile, bwFile)); 		#sort -k1,1N -k2,2n unsrt.bed > srt.bed
		file.remove(bgFile);
		file.remove(tmpFile);
	} else {
		cv = fread(bgFile)
		#cvr = GRanges(seqnames = cv[[1]], ranges = IRanges(start = cv[[2]], end = cv[[3]]), score = cv[[4]])		
		export.bw(cv, bwFile, format = "bw"); #, seqlengths = seqlengths(si[[1]])
	}
	
	invisible(bwFile)
}


gBgFile = "gsub('.bam', paste('.frag', paste(round(fragSizes, 0), collapse = '_'), '.norm', format(normReads, scientific = FALSE), '.bedGraph', sep = ''), bamFiles)";

#Helper function to make a bedgraph file
bamToBedgraph = function(bamFiles., bgFile. = c(gBgFile), normFactors. = c(NA), fragSizes. = c(NA), sbp. = ScanBamParam(scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F)), sigDigits. = c(gSigDigits), message. = c(NA)) {

	#Debug
	#bamFiles. = bamFiles; bgFile = bgFiles[[chr]]; fragSizes. = fragSizes; sigDigits. = gSigDigits; normFactors = normReads / unlist(readCounts); sbp = sbps[[chr]]
	#bamFiles. = bamFiles; bgFile. = bgFiles[[chr]]; fragSizes. = fragSizes; sigDigits. = sigDigits; normFactors. = normReads / unlist(readCounts); sbp. = sbps[[chr]]; message. = chr;
	#bamFiles. = bamFiles; bgFile. = bgFiles[1]; normFactors. = rep(list(normReads / unlist(readCounts)), length(chrs))[1]; fragSizes. = fragSizes; bwFile = bwFile; sbp. = sbps[[1]]; sigDigits. = gSigDigits

	if (!is.na(message.)) print(paste(message., Sys.time()))

	if (bgFile. == gBgFile) bgFile. = eval(parse(text = gBgFile))
	
	#coverage variables
	cvs = RleList(compress = F);
	cv = RleList(compress = F); 

	#Recursively build coverage for each bamFile 
	for (i in 1:length(bamFiles.)) {
	
		#Read in chromosome from each bam file		
		bam = readGAlignments(as.character(bamFiles.[i]), param = sbp.);

		#Extend reads to fragment size
		if (!is.blank(fragSizes.[i])) bam = suppressWarnings(resize(granges(bam), width = fragSizes.[i], fix = "start"));

		#calculate coverage
		cv = list()
		cv[[i]] = coverage(bam)
		for (chr in names(cv[[i]])) if (all(cv[[i]][[chr]]@values == 0)) cv[[i]][[chr]] = NULL

		#normalize coverage to reads
		if (!is.na(normFactors.)) cv[[i]] = round(cv[[i]] * normFactors.[[i]], sigDigits.)

		if (i == 1) cvs = cv[[i]] else cvs = cvs + cv[[i]]

		#clear memory
		rm(bam); gc();
	}

	cvs = cvs / length(bamFiles.); #adjust for # of bamFiles

	#convert coverage Rle object to data table for writing
	dts = NULL
	for (chr in names(cvs)) {
		dt = data.table(chr = chr, start = c(1, cumsum(cvs[[chr]]@lengths)[1:length(cvs[[chr]]@lengths) - 1]), end = cumsum(cvs[[chr]]@lengths), rpm = cvs[[chr]]@values)
		dts = rbind(dts, dt)
	}
	dts = dts[, lapply(.SD, format, scientific = F)]

	#export.bedGraph(cvs, bgFile.)
	fwrite(dts, file = bgFile., sep = "\t", quote = F, row.names = F, col.names = F)
}

writeSiChromSizes = function(si, file, removeChrs = c(NA)) {

	if (!is.na(removeChrs)) inc = !grepl(removeChrs, si@seqnames) else inc = T

	write.table(cbind(si@seqnames[inc], si@seqlengths[inc]), file = file, sep = "\t", quote = F, row.names = F, col.names = F)

}

#### function to get bam pileup values

#### Constants for Hist functions
gHistBin = 10;
gHistOutFile = "paste(bamFiles, '.', gsub('^.*/', '', bedFile), '.range', range, '.bin', bin, '.csv', sep = '')"

makeBamHist = function(bamFiles, bed, outFile = c(gHistOutFile), range, bin = c(gHistBin), fragSizes = c(gFragSize), normReads = c(gNormReads), readCounts = c(NA), sigDigits = c(gSigDigits), flag = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F), threads = c(1), sense = c(any(strand(bed) != "*")), ...) {

	#Debug
	#bamFiles = bamFile; bed = proms; outFile = outFile; range = range; bin = bin; fragSizes = gFragSize; normReads = gNormReads; readCounts = NA; sigDigits = gSigDigits; flag = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F); threads = 24; sense = T
	#bamFiles = bamFile; bed = motifr; outFile = statFile; range = range; bin = bin; fragSizes = fragSizes; normReads = 1e6 / atacSamples$frip.0.01[j]; readCounts = NA; sigDigits = gSigDigits; flag = scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F); threads = threads; sense = T

	#Determine mean width of bed range
	meanWidth = round(mean(width(bed))) - 1

	if (outFile == gHistOutFile)  outFile = eval(parse(text = outFile))

	if (all(is.na(readCounts))) {
		readCounts = list();
		for (i in 1:length(bamFiles)) readCounts[[i]] = sum(idxstatsBam(bamFiles[i])$mapped)
	}

	si = seqinfo(BamFile(bamFiles[1])); 	#get bam genome info
	chrs = as.character(unique(seqnames(bed))); 	#get list of chrs
	chrs = chrs[chrs %in% seqnames(si)]

	#make bed files for each chromosome
	bedChrs = list()
	sbps = list()	
	outFileChrs = paste0(outFile, ".", chrs)
	col.names = c(T, rep(F, length(chrs) - 1)) 
	for (chr in chrs) {
		#build list of beds by chrom
		bedChrs[[chr]] = bed[bed@seqnames == chr];

		#scan bam param for chromosome
		sbps[[chr]] = ScanBamParam(flag = flag, which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))
	}

	hists = mcmapply(bamToHist, bamFiles. = rep(bamFiles, length(bedChrs)), bed. = bedChrs, outFile. = outFileChrs, range. = range, bin. = bin, meanWidth. = meanWidth, fragSizes. = fragSizes, normReads. = normReads, readCounts. = readCounts, sigDigits. = sigDigits, sbp = sbps, chr. = chrs, col.names. = col.names, sense. = sense, mc.cores = threads)

	#concatenate and clean-up
	catFiles(outFileChrs, outFile)
	file.remove(outFileChrs)

	invisible(outFile)
}

#### Helper function to parallelize makebamHist
bamToHist = function(bamFiles., bed., outFile., range., bin. = c(gHistBin), meanWidth. = c(0), fragSizes. = c(gFragSize), normReads. = c(gNormReads), readCounts. = c(NA), sigDigits. = c(gSigDigits), sbp = ScanBamParam(scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = F, isDuplicate = F)), chr. = c(NA), col.names. = c(F), sense. = c(any(strand(bed.) != "*")), ...) {

	#Debug
	#bamFiles. = bamFiles; bed. = bedChrs[[1]]; outFile. = outFileChrs[1]; range. = range; bin. = bin; meanWidth. = meanWidth; fragSizes. = fragSizes; sigDigits. = sigDigits; normReads. = normReads; readCounts. = readCounts; sbp = sbps[[1]]; chr. = chrs[1]; col.names. = col.names[1]; sense. = sense

	print(paste(chr., Sys.time()))

	#binned region to calculate enrichment
	bins = seq(from = -range., to = range. + meanWidth., by = bin.)

	#Create bin labels, normalize distance within bed region 
	binlabs = bins
	binlabs[binlabs > 0 & binlabs <= meanWidth.] = binlabs[binlabs > 0 & binlabs <= meanWidth.] / meanWidth.
	binlabs[binlabs > meanWidth.] = binlabs[binlabs > meanWidth.] - meanWidth. + 1

	if (!is.null(bed.@elementMetadata$name) & length(bed.@elementMetadata$name) == length(unique(bed.@elementMetadata$name))) {
		histRows = bed.@elementMetadata$name
	} else histRows = paste(as.character(bed.@seqnames), start(bed.), end(bed.), strand(bed.), sep = "_");

	#Make hist matrix index
	histIdx = matrix(bins, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);
	histIdx = histIdx + start(bed.)
	histIdx[histIdx <= 0] = NA

	#make hist matrix
	hist = matrix(0, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);

	hists = list();

	#Recursively build coverage for each bamFile 
	for (i in 1:length(bamFiles.)) {

		hists[[i]] = hist;

		if (is.na(fragSizes.[i])) {
			bam = granges(readGAlignments(as.character(bamFiles.[i]), param = sbp));
		} else bam = suppressWarnings(resize(granges(readGAlignments(as.character(bamFiles.[i]), param = sbp)), width = fragSizes.[i], fix = "start"));

		if (is.na(sense.)) {
			#Calculate base level coverage
			cv = coverage(bam)[[chr.]];

			#Calculate the running mean to put into bins
			means = runmean(cv, k = bin.)

			#Get the names of rows for chromosome chr
			chrRows = rownames(hists[[i]])[as.character(seqnames(bed.)) == chr.];

			#Replace regions that are larger than than the genome with NA
			si = seqinfo(BamFile(bamFiles.[i])); 	#get bam genome info	
			hists[[i]][chrRows, ][hists[[i]][chrRows, ] > si@seqlengths[si@seqnames == chr.]] = NA

			#Replace hist rows with mean
			hists[[i]][chrRows, ] = as.numeric(means)[histIdx[chrRows, ]]
	
			rm(cv); rm(means); rm(chrRows); gc()
	
		} else if (!is.na(sense.)) {
			bamPos = bam[strand(bam) == "+"]
			bamNeg = bam[strand(bam) == "-"]

			#Calculate base level coverage
			cvPos = coverage(bamPos)[[chr.]];
			cvNeg = coverage(bamNeg)[[chr.]];

			#Calculate the running mean to put into bins
			meansPos = runmean(cvPos, k = bin.)
			meansNeg = runmean(cvNeg, k = bin.)

			#Get the names of rows for chromosome chr
			chrRowsPos = rownames(hists[[i]])[as.character(seqnames(bed.)) == chr. & as.character(strand(bed.)) == "+"];
			chrRowsNeg = rownames(hists[[i]])[as.character(seqnames(bed.)) == chr. & as.character(strand(bed.)) == "-"];

			#Replace regions that are larger than than the genome with NA
			si = seqinfo(BamFile(bamFiles.[i])); 	#get bam genome info	

			hists[[i]][chrRowsPos, ][hists[[i]][chrRowsPos, ] > si@seqlengths[si@seqnames == chr.]] = NA
			hists[[i]][chrRowsNeg, ][hists[[i]][chrRowsNeg, ] > si@seqlengths[si@seqnames == chr.]] = NA

			#Replace hist rows with mean account for sense
			if (sense.) {
				hists[[i]][chrRowsPos, ] = as.numeric(meansPos)[histIdx[chrRowsPos, ]]
				hists[[i]][chrRowsNeg, ] = as.numeric(meansNeg)[histIdx[chrRowsNeg, ]]
			} else {
				hists[[i]][chrRowsPos, ] = as.numeric(meansNeg)[histIdx[chrRowsPos, ]]
				hists[[i]][chrRowsNeg, ] = as.numeric(meansPos)[histIdx[chrRowsNeg, ]]
			}
									
			rm(bamPos); rm(bamNeg); rm(cvPos); rm(cvNeg); rm(meansPos); rm(meansNeg); rm(chrRowsPos); rm(chrRowsNeg); gc();
		}
		
		rm(bam); gc();
	}

	#add up hists for each bam file
	histSum = NA;
	for (i in 1:length(bamFiles.)) if (all(is.na(histSum))) histSum = hists[[i]] else histSum = histSum + hists[[i]];
		
	#Adjust for strandness
	histSum[as.character(bed.@strand) == "-", ] = histSum[as.character(bed.@strand) == "-", rev(seq_len(ncol(histSum)))];

	colnames(histSum) = binlabs; #label bins correctly	
	#Normalize and round
	histAvg = histSum * normReads. / sum(unlist(readCounts.));
	histAvg = cbind(row = row.names(histAvg), round(histAvg, sigDigits.))
	
	fwrite(as.data.table(histAvg), file = outFile., sep = ",", row.names = F, quote = F, col.names = col.names.) 
	#write.table(histAvg, file = outFile., sep = ",", row.names = F, quote = F, col.names = col.names.) 
	
	rm(histSum); rm(histAvg); gc();
}


#### Constants for Hist functions
gBwHistBin = 10;
gHistOutFile = "paste(bamFiles, '.', gsub('^.*/', '', bedFile), '.range', range, '.bin', bin, '.csv', sep = '')"

makeBwHist = function(bwFiles, bed, outFile = c(gHistOutFile), range, bin = c(gBwHistBin), sigDigits = c(gSigDigits), threads = c(1), sense = c(any(strand(bed) != "*")), ...) {

	#Debug
	#bwFiles = bwFile; bed = motifr; outFile = statFile; range = range; bin = 1; sigDigits = gSigDigits; threads = 5; sense = any(strand(bed) != "*")
		
	#Determine mean width of bed range
	meanWidth = round(mean(width(bed)))

	if (outFile == gHistOutFile)  outFile = eval(parse(text = outFile))

	si = seqinfo(BigWigFile(bwFiles[1])); 	#get bam genome info
	chrs = as.character(unique(seqnames(bed))); 	#get list of chrs
	chrs = chrs[chrs %in% seqnames(si)]

	#make bed files for each chromosome
	bedChrs = list()
	bws = list() # BigWigSelection
	outFileChrs = paste0(outFile, ".", chrs)
	col.names = c(T, rep(F, length(chrs) - 1)) 
	for (chr in chrs) {
		#build list of beds by chrom
		bedChrs[[chr]] = bed[bed@seqnames == chr];

		#bigWigSelections for chromosomes
		#bws[[chr]] = BigWigSelection(GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))
		bws[[chr]] = GRanges(seqnames = chr, ranges = IRanges(start = 1, end = si@seqlengths[si@seqnames == chr]), seqinfo = seqinfo(BigWigFile(bwFiles[1])))
	}

	hists = mcmapply(bamToHist, bwFiles. = rep(bwFiles, bedChrs), bed. = bedChrs, outFile. = outFileChrs, range. = range, bin. = bin, meanWidth. = meanWidth, bws. = bws, sigDigits. = sigDigits, chr. = chrs, col.names. = col.names, sense. = sense, mc.cores = threads)

	#concatenate and clean-up
	catFiles(outFileChrs, outFile)
	file.remove(outFileChrs)

	invisible(outFile)
}

#### Helper function to parallelize makeBwHist
bwToHist = function(bwFiles., bed., outFile., range., bin. = c(gHistBin), meanWidth. = c(0), bws. = BigWigSelection(), sigDigits. = c(gSigDigits), chr. = c(NA), col.names. = c(F), sense. = c(any(strand(bed.) != "*")), ...) {

	#Debug
	#bwFiles. = bwFiles; bed. = bedChrs[[1]]; outFile. = outFileChrs[1]; range. = range; bin. = bin; meanWidth. = meanWidth; bws. = bws[[1]]; sigDigits. = sigDigits; chr. = chrs[1]; col.names. = col.names[1]; sense. = any(strand(bed.) != "*")

	print(paste(chr., Sys.time()))

	#binned region to calculate enrichment
	bins = seq(from = -range., to = range. + meanWidth., by = bin.)

	#Create bin labels, normalize distance within bed region 
	binlabs = bins
	binlabs[binlabs > 0 & binlabs <= meanWidth.] = binlabs[binlabs > 0 & binlabs <= meanWidth.] / meanWidth.
	binlabs[binlabs > meanWidth.] = binlabs[binlabs > meanWidth.] - meanWidth. + 1

	if (!is.null(bed.@elementMetadata$name) & length(bed.@elementMetadata$name) == length(unique(bed.@elementMetadata$name))) {
		histRows = bed.@elementMetadata$name
	} else histRows = paste(as.character(bed.@seqnames), start(bed.), end(bed.), strand(bed.), sep = "_");

	#Make hist matrix index
	histIdx = matrix(bins, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);
	histIdx = histIdx + start(bed.)
	histIdx[histIdx <= 0] = NA

	#make hist matrix
	hist = matrix(NA, ncol = length(bins), nrow = length(histRows), dimnames = list(histRows, bins), byrow = T);

	hists = list();

	#Recursively build coverage for each bamFile 
	for (i in 1:length(bwFiles.)) {

		hists[[i]] = hist;

		bw = import.bw(bwFiles.[i], which = bws.)

		trs = GRanges(seqnames = chr., ranges = IRanges(start = as.vector(t(histIdx)), end = as.vector(t(histIdx))))
		fo = findOverlaps(bw, trs)
		so = subsetByOverlaps(bw, trs)

		hists[[i]][217] = 1
	
		if (is.na(sense.)) {
			#Calculate base level coverage
			cv = coverage(bw)[[chr.]];
			
			#Calculate the running mean to put into bins
			means = runmean(cv, k = bin.)

			#Get the names of rows for chromosome chr
			chrRows = rownames(hists[[i]])[as.character(seqnames(bed.)) == chr.];

			#Replace regions that are larger than than the genome with NA
			si = seqinfo(BamFile(bwFiles.[i])); 	#get bam genome info	
			hists[[i]][chrRows, ][hists[[i]][chrRows, ] > si@seqlengths[si@seqnames == chr.]] = NA

			#Replace hist rows with mean
			hists[[i]][chrRows, ] = as.numeric(means)[histIdx[chrRows, ]]
	
			rm(cv); rm(means); rm(chrRows); gc()
	
		} else if (!is.na(sense.)) {
			bamPos = bam[strand(bam) == "+"]
			bamNeg = bam[strand(bam) == "-"]

			#Calculate base level coverage
			cvPos = coverage(bamPos)[[chr.]];
			cvNeg = coverage(bamNeg)[[chr.]];

			#Calculate the running mean to put into bins
			meansPos = runmean(cvPos, k = bin.)
			meansNeg = runmean(cvNeg, k = bin.)

			#Get the names of rows for chromosome chr
			chrRowsPos = rownames(hists[[i]])[as.character(seqnames(bed.)) == chr. & as.character(strand(bed.)) == "+"];
			chrRowsNeg = rownames(hists[[i]])[as.character(seqnames(bed.)) == chr. & as.character(strand(bed.)) == "-"];

			#Replace regions that are larger than than the genome with NA
			si = seqinfo(BamFile(bwFiles.[i])); 	#get bam genome info	

			hists[[i]][chrRowsPos, ][hists[[i]][chrRowsPos, ] > si@seqlengths[si@seqnames == chr.]] = NA
			hists[[i]][chrRowsNeg, ][hists[[i]][chrRowsNeg, ] > si@seqlengths[si@seqnames == chr.]] = NA

			#Replace hist rows with mean account for sense
			if (sense.) {
				hists[[i]][chrRowsPos, ] = as.numeric(meansPos)[histIdx[chrRowsPos, ]]
				hists[[i]][chrRowsNeg, ] = as.numeric(meansNeg)[histIdx[chrRowsNeg, ]]
			} else {
				hists[[i]][chrRowsPos, ] = as.numeric(meansNeg)[histIdx[chrRowsPos, ]]
				hists[[i]][chrRowsNeg, ] = as.numeric(meansPos)[histIdx[chrRowsNeg, ]]
			}
									
			rm(bamPos); rm(bamNeg); rm(cvPos); rm(cvNeg); rm(meansPos); rm(meansNeg); rm(chrRowsPos); rm(chrRowsNeg); gc();
		}
		
		rm(bam); gc();
	}

	#add up hists for each bam file
	histSum = NA;
	for (i in 1:length(bwFiles.)) if (all(is.na(histSum))) histSum = hists[[i]] else histSum = histSum + hists[[i]];
		
	#Adjust for strandness
	histSum[as.character(bed.@strand) == "-", ] = histSum[as.character(bed.@strand) == "-", rev(seq_len(ncol(histSum)))];

	colnames(histSum) = binlabs; #label bins correctly	
	#Normalize and round
	histAvg = histSum * normReads. / sum(unlist(readCounts.));
	histAvg = cbind(row = row.names(histAvg), round(histAvg, sigDigits.))
	
	fwrite(as.data.table(histAvg), file = outFile., sep = ",", row.names = F, quote = F, col.names = col.names.) 
	#write.table(histAvg, file = outFile., sep = ",", row.names = F, quote = F, col.names = col.names.) 
	
	rm(histSum); rm(histAvg); gc();
}


### annotates a bigwig file for coverage relative to a bed
gListOutFile = "paste(bwFile, '.', gsub('^.*/', '', bedFile), '.range', range, '.csv', sep = '')"
makeBwScoreList = function(bwFile, bed, range = c(NA), outFile = c(gListOutFile), sumFile = c(NA), sigDigits = c(gSigDigits), ...) {

	#Debug
	#bwFile = bwFiles[1]; bed = motifr; range = range; outFile = statFiles[1]; sigDigits = gSigDigits; sumFile = "test.csv"

	if (is.character(bed) && file.exists(bed)) {
		bed = unique(import.bed(bed))
	} else if (!inherits(bed, what = "GRanges")) {
		print("ERROR: bed not bed file or GRanges object");
	}

	if (outFile == gListOutFile) outFile = eval(parse(text = outFile))

	scores = data.frame(chr = NULL, start = NULL, score = NULL, relLoc = NULL)

	si = seqinfo(BigWigFile(bwFile))

	chrs = intersect(si@seqnames, unique(seqnames(bed)))
	
	for (chr in chrs) {

		print(paste("chr", chr, Sys.time()))

		bedChr = bed[bed@seqnames == chr];

		if (!is.na(range)) {
			iRanges = IRanges(start(bedChr) - range, end(bedChr) + range)
			start(iRanges)[start(iRanges) < 1] = 1
			end(iRanges)[end(iRanges) > si@seqlengths[si@seqnames == chr]] = si@seqlengths[si@seqnames == chr]
		} else iRanges = IRanges(1, si@seqlengths[si@seqnames == chr]);

		bw = import.bw(as.character(bwFile), which = GRanges(seqnames = as.character(chr), ranges = iRanges, seqinfo = si))
		
		if (length(bw) > 0) {
			scoresChr = data.frame(chr = chr, start = start(bw), score = round(score(bw), sigDigits), relStart = NA, relEnd = NA, relLoc = NA)

			seqlevels(bw) = chr
			seqlevels(bedChr) = chr;
			dist = distanceToNearest(bw, bedChr)

			bedChrPos = as.character(strand(bedChr[subjectHits(dist)])) == "+"
		
			scoresChr$relStart[bedChrPos] = start(bw)[bedChrPos] - start(bedChr[subjectHits(dist)])[bedChrPos];
			scoresChr$relStart[!bedChrPos] = end(bedChr[subjectHits(dist)])[!bedChrPos] - start(bw)[!bedChrPos];

			scoresChr$relEnd[bedChrPos] = start(bw)[bedChrPos] - end(bedChr[subjectHits(dist)])[bedChrPos];
			scoresChr$relEnd[!bedChrPos] = start(bedChr[subjectHits(dist)])[!bedChrPos] - start(bw)[!bedChrPos];	

			#set relLoc (relative location) to distance. Scale bw elements with bed files from 0 to 1
			beforeBed = scoresChr$relStart <= 0
			scoresChr$relLoc[beforeBed] = scoresChr$relStart[beforeBed];

			#if bw location is within the bed file then scale the relLoc from 0 (start) to 1 (end) directionally
			inBed = scoresChr$relStart > 0 & scoresChr$relEnd < 0
			scoresChr$relLoc[inBed] = scoresChr$relStart[inBed] / (scoresChr$relStart[inBed] - scoresChr$relEnd[inBed]);

			afterBed = scoresChr$relEnd >= 0
			scoresChr$relLoc[afterBed] = scoresChr$relEnd[afterBed];

			scores = rbind(scores, scoresChr)

			rm(bw); rm(scoresChr); gc();
		}
	}

	fwrite(scores, outFile, sep = "\t", row.names = F)


	if (!is.na(sumFile)) {

		bedSize = mean(width(bed))
		#expand for motif size
		scores$relLoc[scores$relLoc >= 1] = scores$relLoc[scores$relLoc >= 1] + bedSize
		scores$relLoc[scores$relLoc >= 0 & scores$relLoc < 1] = scores$relLoc[scores$relLoc >= 0 & scores$relLoc < 1] * bedSize

		mean = aggregate(scores$score, by = list(scores$relLoc), FUN = mean, na.rm = T)
		setnames(mean, c("pos", "mean"))

		ct = aggregate(scores$score, by = list(scores$relLoc), FUN = length)
		setnames(ct, c("pos", "ct"))

		mean = merge(mean, ct, by = "pos")
		mean$mean = format(mean$mean, digits = sigDigits)

		fwrite(mean, sumFile, sep = ",", row.names = F)
		rm(mean); rm(stat); rm(ct); gc()
	}

	invisible(scores)
}


annotBedtoTxDbGene = function(bed, tx, org = c(NULL), prefix = c(NULL), promUp = 2500, promDown = 0, org.keytype = c("ENTREZID"), ...) {

	#Debug
	#bed = tabr; tx = tx; org = NULL; prefix = "GRCh37"; promUp = 2500; promDown = 0; org.keytype = "ENTREZID";
	#bed = peaksBed; tx = tx; org = org; prefix = "GRCh37"; promUP = 2500; promDown = 0;
	#bed = trLocs; tx = tx; org = NULL; prefix = "GRCh37"; promUp = 2500; promDown = 0; org.keytype = "ENTREZID"

	#bed = rnar[95254]; tx = tx; org = NULL; prefix = "GRCh37"; promUp = 2500; promDown = 0; 


	#Columns to add to bed Granges
	annotCols = paste(prefix, c("gene", "promoter", "exon", "intron", "5utr", "3utr", "ts", "tsKg", "tsDist", "tsStart", "tsEnd", "tsStrand", "tsRelPos", "tss", "tssKg", "tssDist", "tssRelPos", "tssStart", "tssEnd", "tssStrand"), sep = ".");
	bed@elementMetadata[annotCols] = rep(NA, dim(bed@elementMetadata)[1])

	#get Ts database
	ts = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")); #, ...); #get transcripts

	#Determine the the closest transcript
	closeTs = distanceToNearest(bed, ts, ignore.strand = T)

	#Assign ts name to annotation column (AnnotCol) "tx"
	bed@elementMetadata[[paste0(prefix, ".ts")]][queryHits(closeTs)] = ts@elementMetadata$TXNAME[subjectHits(closeTs)]
 	bed@elementMetadata[paste0(prefix, ".tsKg")] = select(tx, keytype = "TXNAME", keys = as.character(bed@elementMetadata[[paste0(prefix, ".ts")]]), columns = "GENEID")$GENEID

	bed@elementMetadata[[paste0(prefix, ".tsStart")]][queryHits(closeTs)] = start(ts)[subjectHits(closeTs)]
	bed@elementMetadata[[paste0(prefix, ".tsEnd")]][queryHits(closeTs)] = end(ts)[subjectHits(closeTs)]
	bed@elementMetadata[[paste0(prefix, ".tsStrand")]][queryHits(closeTs)] = as.character(strand(ts)[subjectHits(closeTs)])

	 #Determine the orientation
	sign = sign(start(bed[queryHits(closeTs)]) - start(ts[subjectHits(closeTs)]))
	sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"] = -sign[as.character(strand(ts[subjectHits(closeTs)])) == "-"]; #adjust for strand
	
	#Assign ts dist	
	bed@elementMetadata[[paste(prefix, "tsDist", sep = ".")]][queryHits(closeTs)] = closeTs@elementMetadata$distance * sign

	#Determine ts rel pos
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)] = closeTs@elementMetadata$distance * sign

	#Determine if the bed is within or after a transcript, if so scale it to the transcript
	inTx = bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)] == 0
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)][inTx] = (start(bed[queryHits(closeTs)][inTx]) - start(ts[subjectHits(closeTs)])[inTx]) / (end(ts[subjectHits(closeTs)])[inTx] - start(ts[subjectHits(closeTs)])[inTx])
	bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)][inTx & as.character(strand(ts[subjectHits(closeTs)])) == "-"] = (1 - bed@elementMetadata[[paste0(prefix, ".tsRelPos")]][queryHits(closeTs)][inTx & as.character(strand(ts[subjectHits(closeTs)])) == "-"])

	tss = transcripts(tx, columns = c("TXID", "TXNAME", "GENEID")); #, ...);	
	#make strand-specific ts range object
	end(tss)[as.character(strand(tss)) == "+"] = start(tss)[as.character(strand(tss)) == "+"]
	start(tss)[as.character(strand(tss)) == "-"] = end(tss)[as.character(strand(tss)) == "-"]

	#Determine the relative location
	closeTss = distanceToNearest(bed, tss, ignore.strand = T)
	#closeTss = closeTss[!is.na(subjectHits(closeTss)), ]
	
	#Assign ts name to annotation column (AnnotCol) "ts"
	bed@elementMetadata[[paste0(prefix, ".tss")]][queryHits(closeTs)] = tss@elementMetadata$TXNAME[subjectHits(closeTss)]
	bed@elementMetadata[[paste0(prefix, ".tssKg")]] = select(tx, keytype = "TXNAME", keys = bed@elementMetadata[[paste0(prefix, ".tss")]], columns = "GENEID")$GENEID

	bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)] = start(ts)[subjectHits(closeTss)]
	bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)] = end(ts)[subjectHits(closeTss)]
	bed@elementMetadata[[paste0(prefix, ".tssStrand")]][queryHits(closeTss)] = as.character(strand(ts)[subjectHits(closeTss)])

	#Determine the orientation
	sign = sign(start(bed[queryHits(closeTss)]) - start(tss[subjectHits(closeTss)]))
	sign[as.character(strand(tss[subjectHits(closeTss)])) == "-"] = -sign[as.character(strand(tss[subjectHits(closeTss)])) == "-"]; #adjust for strand

	#Assign tss dist	
	bed@elementMetadata[[paste0(prefix, ".tssDist")]][queryHits(closeTss)] = closeTss@elementMetadata$distance * sign

	#Determine the position relative to the closest tss
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] = closeTss@elementMetadata$distance * sign

	#Determine if the bed is within or after a transcript, if so scale it to the transcript
	inTx = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] >= 0 & bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] < (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)])
	afterTx = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)] >= (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)])
	
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][inTx] = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][inTx] / (bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)][inTx] - bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)][inTx])
	bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][afterTx] = bed@elementMetadata[[paste0(prefix, ".tssRelPos")]][queryHits(closeTss)][afterTx] - bed@elementMetadata[[paste0(prefix, ".tssEnd")]][queryHits(closeTss)][afterTx] + bed@elementMetadata[[paste0(prefix, ".tssStart")]][queryHits(closeTss)][afterTx]

	#Assign gene overlap
	genes = suppressWarnings(genes(tx, columns = c("GENEID")))
	overlaps = findOverlaps(bed, genes)
	overlaps = cbind(as.data.frame(overlaps), id = as.character(genes@elementMetadata$GENEID[subjectHits(overlaps)]))
	overlaps = data.table(unique(overlaps[, c(1, 3)]))
	overlaps = overlaps[!is.na(overlaps$id), ]
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "gene", sep = ".")]][ag$queryHits] = ag$id

	#Assign promoter overlap
	proms = suppressWarnings(promoters(tx, upstream = promUp, downstream = promDown, columns = c("GENEID")))
	overlaps = findOverlaps(bed, proms)
	overlaps = cbind(as.data.frame(overlaps), id = as.character(proms@elementMetadata$GENEID[subjectHits(overlaps)]))
	overlaps = data.table(unique(overlaps[, c(1, 3)]))
	overlaps = overlaps[!is.na(overlaps$id), ]
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "promoter", sep = ".")]][ag$queryHits] = ag$id

	#Assign exon overlap
	exons = exonsBy(tx, by = "tx")
	if (exists("vals") && all(!is.na(vals$tx_name))) exons = exons[names(exons) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset Exon")
	overlaps = findOverlaps(bed, exons)
	overlaps = cbind(as.data.frame(overlaps), id = names(exons)[subjectHits(overlaps)]); 
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXID)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[, c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "exon", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign intron overlap
	introns = intronsByTranscript(tx, use.name = T)
	if (exists("vals") && all(!is.na(vals$tx_name))) introns = introns[names(introns) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset Intron")
	overlaps = findOverlaps(bed, introns)
	overlaps = cbind(as.data.frame(overlaps), id = names(introns)[subjectHits(overlaps)]); 
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "intron", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign 5'utr classification
	fiveUtr = fiveUTRsByTranscript(tx, use.name = T)
	if (exists("vals") && all(!is.na(vals$tx_name))) fiveUtr = fiveUtr[names(fiveUtr) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset 5' UTR")
	overlaps = findOverlaps(bed, fiveUtr)
	overlaps = cbind(as.data.frame(overlaps), id = names(fiveUtr)[subjectHits(overlaps)]); 
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "5utr", sep = ".")]][ag$queryHits] = as.character(ag$id)

	#Assign 3'utr classification
	threeUtr = threeUTRsByTranscript(tx, use.name = T)
	if (exists("vals") && all(!is.na(vals$tx_name))) threeUtr = threeUtr[names(threeUtr) %in% vals$tx_name] else if (exists("vals")) stop("Need vals$tx_name to subset 3' UTR")
	overlaps = findOverlaps(bed, threeUtr)
	overlaps = cbind(as.data.frame(overlaps), id = names(threeUtr)[subjectHits(overlaps)]); 
	overlaps$id = as.character(ts@elementMetadata$GENEID[match(overlaps$id, ts@elementMetadata$TXNAME)]); #switch transcript id with gene id
	overlaps = data.table(unique(overlaps[!is.na(overlaps$id), c(1, 3)])); #Remove duplicates
	ag = overlaps[, list(id = paste(id, collapse = ";")), by = c('queryHits')]
	bed@elementMetadata[[paste(prefix, "3utr", sep = ".")]][ag$queryHits] = as.character(ag$id)

	if (!is.null(org)) {
		bed@elementMetadata[paste0(prefix, ".tsSymbol")] = rep(NA, length(bed))
		sym = unique(select(org, keytype = org.keytype, keys = as.character(bed@elementMetadata[[paste0(prefix, ".tsKg")]]), columns = c(org.keytype, "SYMBOL")))
		sym = aggregate(sym$SYMBOL, by = list(sym[[org.keytype]]), FUN = paste, collapse = "|")

		bed@elementMetadata[[paste0(prefix, ".tsSymbol")]] = sym$x[match(bed@elementMetadata[[paste0(prefix, ".tsKg")]], sym$Group.1)]
			
		bed@elementMetadata[paste0(prefix, ".tssSymbol")] = rep(NA, length(bed))
		sym = unique(select(org, keytype = org.keytype, keys = as.character(bed@elementMetadata[[paste0(prefix, ".tssKg")]]), columns = c(org.keytype, "SYMBOL")))
		sym = aggregate(sym$SYMBOL, by = list(sym[[org.keytype]]), FUN = paste, collapse = "|")
		bed@elementMetadata[[paste0(prefix, ".tssSymbol")]] = sym$x[match(bed@elementMetadata[[paste0(prefix, ".tssKg")]], sym$Group.1)]

	}		
 
	return(bed);
}

annotBedtoGffStart = function(bed, gff, prefix = c(NULL)) {

	#Debug
	#bed = bed; gff = gff; prefix = "test"; fileName = NA;

	#Columns to add to bed Granges
	gffCols = paste(prefix, c("start", "end", "strand", "Name", "distStart", "distEnd"), sep = ".");

	bed@elementMetadata[gffCols] = rep(NA, dim(bed@elementMetadata)[1])

	# get the chromosomes of the bed bed file
	chrs = as.character(seqnames(bed)@values);

	for (chr in chrs) {

		print(paste(chr, Sys.time()))
		
		closeGffs = NULL;
		
		#limit data sets to relevant chromsome chr	
		bedChr = bed[bed@seqnames == chr];
		gffChr = gff[gff@seqnames == chr];

		#set up strand specific start / ends
		gffStarts = start(gffChr);
		gffEnds = end(gffChr);
		gffStarts[as.character(strand(gffChr)) == "-"] = end(gffChr)[as.character(strand(gffChr)) == "-"];
		gffEnds[as.character(strand(gffChr)) == "-"] = start(gffChr)[as.character(strand(gffChr)) == "-"];

		#set up bed starts / ends
		bedStarts = start(bedChr);
		bedStarts[as.character(strand(bedChr)) == "-"] = end(bedChr)[as.character(strand(bedChr)) == "-"];
		bedEnds = end(bedChr);
		bedEnds[as.character(strand(bedChr)) == "-"] = start(bedChr)[as.character(strand(bedChr)) == "-"];

		centers = (start(bedChr) + end(bedChr)) / 2;

		for (i in 1:length(bedChr)) closeGffs = c(closeGffs, which(abs(centers[i] - gffStarts) == min(abs(centers[i] - gffStarts)))[1]);

		bed@elementMetadata[bed@seqnames == chr, gffCols[1]] = gffStarts[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[2]] = gffEnds[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[3]] = as.character(strand(gffChr[closeGffs]))
		bed@elementMetadata[bed@seqnames == chr, gffCols[4]] = gffChr@elementMetadata$Name[closeGffs];
		bed@elementMetadata[bed@seqnames == chr, gffCols[5]] = (bedStarts - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))
		bed@elementMetadata[bed@seqnames == chr, gffCols[6]] = (bedEnds - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))

	}

	return(bed);
}

annotBedtoGffEnd = function(bed, gff, prefix = c(NULL)) {

	#Debug
	#bed = bed; gff = gff; prefix = "test"; fileName = NA;

	#Columns to add to bed Granges
	gffCols = paste(prefix, c("start", "end", "strand", "Name", "distStart", "distEnd"), sep = ".");

	bed@elementMetadata[gffCols] = rep(NA, dim(bed@elementMetadata)[1])

	# get the chromosomes of the bed bed file
	chrs = as.character(seqnames(bed)@values);

	for (chr in chrs) {

		print(paste(chr, Sys.time()))
		
		closeGffs = NULL;
		
		#limit data sets to relevant chromsome chr	
		bedChr = bed[bed@seqnames == chr];
		gffChr = gff[gff@seqnames == chr];

		#set up strand specific start / ends
		gffStarts = start(gffChr);
		gffEnds = end(gffChr);
		gffStarts[as.character(strand(gffChr)) == "-"] = end(gffChr)[as.character(strand(gffChr)) == "-"];
		gffEnds[as.character(strand(gffChr)) == "-"] = start(gffChr)[as.character(strand(gffChr)) == "-"];

		#set up bed starts / ends
		bedStarts = start(bedChr);
		bedStarts[as.character(strand(bedChr)) == "-"] = end(bedChr)[as.character(strand(bedChr)) == "-"];
		bedEnds = end(bedChr);
		bedEnds[as.character(strand(bedChr)) == "-"] = start(bedChr)[as.character(strand(bedChr)) == "-"];

		centers = (start(bedChr) + end(bedChr)) / 2;

		for (i in 1:length(bedChr)) closeGffs = c(closeGffs, which(abs(centers[i] - gffEnds) == min(abs(centers[i] - gffEnds)))[1]);

		bed@elementMetadata[bed@seqnames == chr, gffCols[1]] = gffStarts[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[2]] = gffEnds[closeGffs]
		bed@elementMetadata[bed@seqnames == chr, gffCols[3]] = as.character(strand(gffChr[closeGffs]))
		bed@elementMetadata[bed@seqnames == chr, gffCols[4]] = gffChr@elementMetadata$Name[closeGffs];
		bed@elementMetadata[bed@seqnames == chr, gffCols[5]] = (bedStarts - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))
		bed@elementMetadata[bed@seqnames == chr, gffCols[6]] = (bedEnds - gffStarts[closeGffs]) * unlist(lapply(as.character(strand(gffChr)[closeGffs]) == "-", function(x) (if (x) -1 else 1)))

	}

	return(bed);
}

#Function to annotate a bed object to another bed object 
annotBedtoBed = function(bedA, bed, prefix = c(NULL)) {

	#Debug
	#bedA = covR; bed = bedAnnot; prefix = "cpgi";
	
	chrs = as.character(intersect(seqnames(bed), seqnames(bedA)));
	bedACom = bedA[seqnames(bedA) %in% chrs]
	bedCom = bed[seqnames(bed) %in% chrs]

	seqlevels(bedACom) = chrs;
	seqlevels(bedCom) = chrs;

	#Columns to add to bed Granges
	bedCols = paste(prefix, c("start", "end", "strand", "name", "dist", "relDist"), sep = ".");
	bedACom@elementMetadata[bedCols] = rep(NA, dim(bedACom@elementMetadata)[1])

	closeBed = distanceToNearest(bedACom, bedCom)

	#set the closest 
	bedACom@elementMetadata[, bedCols[1]] = start(bedCom)[closeBed@subjectHits]
	bedACom@elementMetadata[, bedCols[2]] = end(bedCom)[closeBed@subjectHits]
	bedACom@elementMetadata[, bedCols[3]] = as.character(strand(bedCom)[closeBed@subjectHits])

	if (!is.null(bedCom@elementMetadata$name)) bedACom@elementMetadata[, bedCols[4]] = bedCom@elementMetadata$name[closeBed@subjectHits];

	bedACom@elementMetadata[, bedCols[5]] = closeBed@elementMetadata$distance

	#orient distance relative to location and strand (e.g. upstream or downstream)
	upstream = (bedACom@elementMetadata[, bedCols[1]] > end(bedACom) & (bedACom@elementMetadata[bedCols[3]][[1]] == "*" | bedACom@elementMetadata[bedCols[3]][[1]] == "+")) | (bedACom@elementMetadata[, bedCols[2]] < start(bedACom) & bedACom@elementMetadata[bedCols[3]][[1]] == "-") 
	bedACom@elementMetadata[upstream, bedCols[5]] = -bedACom@elementMetadata[upstream, bedCols[5]]

	#Set the relative distance 
	bedACom@elementMetadata[, bedCols[6]] = bedACom@elementMetadata[, bedCols[5]]
	overlapPos = bedACom@elementMetadata[, bedCols[6]] == 0 & (bedACom@elementMetadata[, bedCols[3]] == "+" | bedACom@elementMetadata[bedCols[3]][[1]] == "*")
	overlapNeg = bedACom@elementMetadata[, bedCols[6]] == 0 & bedACom@elementMetadata[, bedCols[3]] == "-"
	bedACom@elementMetadata[overlapPos, bedCols[6]] = (start(bedACom)[overlapPos] - bedACom@elementMetadata[overlapPos, bedCols[1]]) / (bedACom@elementMetadata[overlapPos, bedCols[2]] - bedACom@elementMetadata[overlapPos, bedCols[1]])
	bedACom@elementMetadata[overlapNeg, bedCols[6]] = (bedACom@elementMetadata[overlapNeg, bedCols[2]] - start(bedACom)[overlapNeg]) / (bedACom@elementMetadata[overlapNeg, bedCols[2]] - bedACom@elementMetadata[overlapNeg, bedCols[1]])
	
	bedA@elementMetadata[bedCols] = rep(NA, dim(bedA@elementMetadata)[1])
	for (col in bedCols) bedA@elementMetadata[seqnames(bedA) %in% chrs, col] = bedACom@elementMetadata[col][[1]]

	return(bedA);
}

gMetadata = c("type", "Name");

bedOverlapGff = function(bed, gff, prefix = c(NULL), metadata = gMetadata) {

	#Debug
	#bed = bed; gff = gff; prefix = "test"; metadata = gMetadata;

	#Columns to add to bed Granges
	bedCols = paste(prefix, metadata, sep = ".");
	bed@elementMetadata[bedCols] = rep(NA, dim(bed@elementMetadata)[1])

	laps = findOverlaps(bed, gff)

	for (i in unique(laps@queryHits)) for (j in 1:length(bedCols)) bed@elementMetadata[i, bedCols[j]] = paste(gff@elementMetadata[metadata[j]][[1]][laps@subjectHits[laps@queryHits == i]], collapse = ";")

	return(bed);
}


#Function to consolidate bed intersection
bedInt = function(bed1, bed2, outFile = c(NA)) {

	if (class(bed1) == "character") bed1 = import.bed(bed1)
	if (class(bed2) == "character") bed2 = import.bed(bed2)

	int = intersect(bed1, bed2)

	score1 = aggregate(bed1@elementMetadata$score[findOverlaps(bed1, int)@queryHits], by = list(findOverlaps(bed1, int)@subjectHits), FUN = mean, na.rm = T)$x
	score2 = aggregate(bed2@elementMetadata$score[findOverlaps(bed2, int)@queryHits], by = list(findOverlaps(bed2, int)@subjectHits), FUN = mean, na.rm = T)$x

	score = apply(cbind(score1, score2), MARGIN = 1, mean, na.rm = TRUE)
	values(int) = data.frame(score = score);

	if (!is.na(outFile)) write.table(int, outFile = outFile, sep = "\t", row.names = F);

	return(int);
}



#Function to consolidate bed union
bedUnion = function(bed1File, bed2File) {

	#Debug

	bed1 = import.bed(bed1File)
	bed2 = import.bed(bed2File)

	bed = union(bed1, bed2)

	score1 = rep(NA, length(bed));
	score2 = rep(NA, length(bed));

	score1[findOverlaps(bed1, bed)@subjectHits] = bed1@elementMetadata$score
	score2[findOverlaps(bed2, bed)@subjectHits] = bed2@elementMetadata$score
	
	score = apply(cbind(score1, score2), MARGIN = 1, mean, na.rm = TRUE)

	name1 = rep(NA, length(bed));
	name2 = rep(NA, length(bed));

	name1[findOverlaps(bed1, bed)@subjectHits] = bed1@elementMetadata$name
	name2[findOverlaps(bed2, bed)@subjectHits] = bed2@elementMetadata$name

	name = paste(name1, name2, sep = "_")

	values(bed) = data.frame(score = score, name = name);

	return(bed);
}


#Function to consolidate bed union
exportBed = function(bed, con) {

	bed = amp; con = "test.bed"

	#Debug
	out = as.data.frame(bed);
	if (any(grepl("strand", names(out)))) out$strand = gsub("\\*", "\\.", out$strand)
	write.table(out, con, sep = "\t", quote = F, col.names = F, row.names = F)	
}

#function to shuffle regions in a given genome
#regions = GRanges object
#genome = BSgenome object
shuffleReg = function(regions, genome, resolution = c(1e10)) {

	chr = sample(seqlevels(genome), length(regions), replace = T, prob = seqlengths(genome)) 
	s = sample.int(resolution, length(regions)) / resolution
	start = round(s * seqlengths(genome)[chr] - width(regions), 0)
	random = GRanges(seqnames = chr, ranges = IRanges(start = start, end = start + width(regions)))

	return(random)
}

 
