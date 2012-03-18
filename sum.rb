#!/usr/bin/ruby

T = 2000
Binsize = 10
DT = 0.001

N = 102400
Input = "gr.spk.a1"
Output = "sum.dat"

ary = Array.new(T/Binsize)

IO.foreach(Input){|l|
  t, n = l.chomp.split
  ary[t.to_i/Binsize] = ary[t.to_i/Binsize].to_i + 1
}

open(Output, "w"){|o|
  ary.each_with_index{|n, i|
    o.puts "#{DT*i.to_f*Binsize} #{n.to_f/(N.to_f)}"
  }
}
