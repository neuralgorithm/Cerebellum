#!/usr/bin/ruby

initial = true
initial_modulation = 1.0
open("gainchange.dat", "w"){|o|
  (0..300).step(10){|i|
    system "./psth.rb #{i} 16"
    fn = "#{i}_16.dat"
    system "cp #{fn} tmp; gnuplot fit.gp > tmp.dat 2>&1"
    IO.foreach("tmp.dat"){|l|
      if l =~ /\A\s+a = /
        modulation = l.gsub(/\A\s+a = /, "").gsub(/\s+\Z/, "").to_f.abs
        if initial
          initial_modulation = modulation
          initial = false
        end
        o.puts "#{i} #{modulation/initial_modulation}"
      end
    }
  }
}
system "rm -f tmp tmp.dat"
