#! /bin/env ruby

require 'getoptlong'

exp_file,prediction_file=nil

#########################################################
class MICRORNA_target_interact_list
  attr_accessor :file
  attr_reader :interaction_list
  def initialize(file)
    @file=file
    @interaction_list=Hash.new
  end
  def read_info
    fh=File.open(@file)
    while(line=fh.gets) do
      line.chomp!
      target, *microRNAs = line.split("\t")
      microRNAs.each do |microRNA|
        microRNA.upcase!
        @interaction_list[[target,microRNA].join('-')]=1
      end
    end
  end
end

#########################################################
opts = GetoptLong.new(
  ['--exp_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--prediction_file', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--exp_file'
      exp_file=value
    when '--prediction_file'
      prediction_file=value
  end
end

#########################################################
exp_obj=MICRORNA_target_interact_list.new(exp_file)
exp_obj.read_info
prediction_obj=MICRORNA_target_interact_list.new(prediction_file)
prediction_obj.read_info

#exp_obj.interaction_list.each_key do |interaction|
#  puts interaction
#end

#intersect_array = prediction_obj.interaction_list.keys - exp_obj.interaction_list.keys

intersect_array = prediction_obj.interaction_list.keys & exp_obj.interaction_list.keys

puts intersect_array.size, exp_obj.interaction_list.keys.size
puts (intersect_array.size.to_f/exp_obj.interaction_list.keys.size)

