#!/usr/bin/env ruby
#Merge the two filtered gene lists according to common genes.

ifile1 = ARGV[0] # scripts/Mus_musculus.NCBIM37.63.gtf
ifile2 = ARGV[1] # Brain/DE.w.limma.voom/Sleep.vs.SleepDep.w_LimmaVoom_DEG_results.txt 

genelibrary = {}
File.open(ifile1,"r").each_line do |line|
	line.strip!
	genetype = line.split("\t")[1]
	geneid = line.split("\t")[8].split(";")[0].split("\"")[1]
	genelibrary[geneid] = genetype
end

File.open(ifile2,"r").each_line do |line|
	line.strip!
	if line.start_with?("\"id\"")
		puts [line, "genetype"].join("\t")
		next
	end
	geneid = line.split("\t")[0].split("\"")[1] # 55 for intron, 0 for gene
	if geneid.include?(",")
		geneids = geneid.split(",")
		geneids.each do |id|
			if genelibrary.has_key?(id)
				puts [line, genelibrary[id]].join("\t")
			end
		end
	else
		if genelibrary.has_key?(geneid)
			puts [line, genelibrary[geneid]].join("\t")
		end
	end
end

