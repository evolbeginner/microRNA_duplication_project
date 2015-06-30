#! /bin/env ruby

require 'bio'
require 'getoptlong'

##############################################################
def get_corename_of_miRNA(miRNA_name)
  return_value = false
  if miRNA_name =~ /^(?:[^-]+) [-] (mir\d+)(\D*)/ix then # miRNA_name e.g. : ath-MIR156e
    return_value=$1
  end
  return(return_value)
end

##############################################################
mature_miRNAs_file=File.expand_path("~/project/microRNA/sequence/miRNA_seq/mature.fa")
search_4_result_dir=''
present={}
excluded=[]
included=[]
out_dir='.'
prefix=''
miRNA_diver_result_file_arr=[]
is_no_output_list=FALSE

opts = GetoptLong.new(
  [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
  [ '--search_4_result_dir', GetoptLong::REQUIRED_ARGUMENT],
  [ '--excluded', GetoptLong::REQUIRED_ARGUMENT],
  [ '--miRNA_diver_result_file', GetoptLong::REQUIRED_ARGUMENT],
  [ '--out_dir', GetoptLong::REQUIRED_ARGUMENT],
  [ '--prefix', GetoptLong::REQUIRED_ARGUMENT],
  [ '--no_output_list', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, arg|
  case opt
    when '--search_4_result_dir'
      search_4_result_dir=arg
    when '--excluded'
      excluded.push(arg)
    when '--prefix'
      prefix=arg
    when '--out_dir'
      out_dir=arg
    when '--miRNA_diver_result_file'
      miRNA_diver_result_file_arr.push(arg)
    when '--no_output_list'
      is_no_output_list=1
  end
end

File.directory?(out_dir) ? () : (Dir.mkdir(out_dir))

###############################################################
present["duplicates"]={}
present["singletons"]={}
all_gene={}

Dir.foreach(search_4_result_dir){|dir1| # dir1, e.g. : A_lyrata
  dir1 =~ /^\./ and next
  FileTest::exists?(dir1+'/'+"duplicates.list") and next
  file_names_hash={}
  types_arr=["duplicates", "singletons"]
  types_arr.each{|x|
    present[x][dir1]={}
    file_name=File.join(search_4_result_dir, dir1, x+".list")
    # file_names_hash[x] = file_name
    fh = File.open(file_name)
    while line = fh.gets do
      line=line.chomp
      line=get_corename_of_miRNA(line)
      next if ! line
      present[x][dir1].merge!({line=>1})
      all_gene[line]=1
    end
  }
}

included_duplicates = present["duplicates"].keys.select{|key| key if ! excluded.include?(key) }
present["singletons"].delete_if{|key,value| excluded.include?(key) }
present["duplicates"].delete_if{|key,value| excluded.include?(key) }

included_miRNAs={}
included.each do |species_name|
  present["duplicates"][species_name].each do |miRNA_name|
    included_miRNAs[miRNA_name]=1
  end
end

############################################################
num_of_gene={}
num_of_gene["singletons"]={}; num_of_gene["duplicates"]={}
present.each_key do |gene_type| # e.g. duplicates or singletons
  present[gene_type].each_key do |species_name|
    present[gene_type][species_name].each_key do |miRNA_name|
      if gene_type == 'singletons' then
        if ! present["duplicates"][species_name].include?(miRNA_name) then
          (num_of_gene[gene_type][miRNA_name].class==Fixnum) ? (num_of_gene[gene_type][miRNA_name]+=1) : (num_of_gene[gene_type][miRNA_name]=1)
        end
      else
          (num_of_gene[gene_type][miRNA_name].class==Fixnum) ? (num_of_gene[gene_type][miRNA_name]+=1) : (num_of_gene[gene_type][miRNA_name]=1)
      end
    end
  end
end

final_list_of_singletons = num_of_gene["singletons"].select {|x,y| y == present["singletons"].size}

######  get all (plant) miRNAs names and eliminate those from the final list of singletons ######
fh = Bio::FlatFile.open(mature_miRNAs_file)
fh.each_entry do |f|
  if f.definition =~ /([a-z]{3})\-miR(\d+)/ then
    orgn_name=$1
    miRNA_num=$2
    if orgn_name !~ /^ath$|^aly$/ then
      miRNA_full_name="MIR"+miRNA_num
      final_list_of_singletons.delete(miRNA_full_name)
    end
  end
end

###### Outputting ######
if ! is_no_output_list then
  sin_fh = File.open([out_dir + '/' + prefix, "singletons_seq.list"].join('-'), 'w')
  dup_fh = File.open([out_dir + '/' + prefix, "duplicates_seq.list"].join('-'), 'w')
  final_list_of_singletons.each {|x| sin_fh.print x[0],"\n"}
  all_gene.each_key{|gene_name| dup_fh.puts gene_name if ! final_list_of_singletons.include?(gene_name)}
end


# num_of_gene["duplicates"].each_key{|miRNA_name| puts miRNA_name}
##################################################################
miRNA_diver_result_file_arr.each do |file|
  fh = File.open(file, 'r')
  while line=fh.gets do
    a=FALSE
    line=line.chomp
    line_arr = line.split("\t")
    line_arr.values_at(1..line_arr.size-1).each do |gene_name|
      gene_name.upcase!
      gene_name = gene_name.sub(/\D?$/, "")
      if final_list_of_singletons.include?(gene_name) then
        print gene_name + "\t"
        a=TRUE
      end
    puts if a
    end
  end  
end

