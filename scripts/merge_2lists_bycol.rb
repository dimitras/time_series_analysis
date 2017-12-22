#!/usr/bin/env ruby
#Merge tables according to column.

#ruby merge_2lists_bycol.rb voila.SS_3hr.vs.SD_3hr/SS_3hr_SD_3hr.deltapsi_deltapsi.tsv DE.w.limma.voom/Sleep.vs.SleepDep.3h.w_LimmaVoom_DEG_results.txt > SS_3hr_SD_3hr_spliced_in_degs.txt
#ruby merge_2lists_bycol.rb voila.SD_3hr.vs.SD_12hr/SD_3hr_SD_12hr.deltapsi_deltapsi.tsv DE.w.limma.voom/Sleep.vs.SleepDep.w_LimmaVoom_DEG_results.txt > SD_3hr_SD_12hr_spliced_in_degs.txt

#ruby merge_2lists_bycol.rb voila.SS_3hr.vs.SD_3hr/SS_3hr_SD_3hr.deltapsi_deltapsi.tsv master_list_of_gene_counts_MIN.sense.Brain.txt > SS_3hr_SD_3hr_spliced_in_fulllist.txt
#ruby merge_2lists_bycol.rb voila.SD_3hr.vs.SD_12hr/SD_3hr_SD_12hr.deltapsi_deltapsi.tsv master_list_of_gene_counts_MIN.sense.Brain.txt > SD_3hr_SD_12hr_spliced_in_fulllist.txt

ifile1 = ARGV[0]
ifile2 = ARGV[1]

list1 = {}
File.open(ifile1,"r").each_line do |line|
        line.strip!
        if line.start_with?("#Gene Name")
                next
        end
        id = line.split("\t")[1]
        list1[id] = "0.1"
end

File.open(ifile2,"r").each_line do |line|
        line.strip! 
        if line.start_with?("id")
                puts [line.split("\t")[0..54], "spliced"].join("\t") #57 for degs list, 54 for full port list
                next
        end
        id = line.split("\t")[0]
        spliced = nil
        if !list1.has_key?(id)
                spliced = 1
        elsif list1.has_key?(id)
                spliced = 0.1           
        end
        puts [line.split("\t")[0..54], spliced].join("\t") #57 for degs list, 54 for full port list
end