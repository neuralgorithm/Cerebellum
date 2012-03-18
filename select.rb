#!/usr/bin/ruby

T = 2000
Binsize = 10
DT = 0.001
N = 100

indices = [136, 18, 86, 117, 94]
files = Array.new(indices.size)

indices.each_with_index{|n, i|
  files[i] = open("gr.#{n}", "w")
}

(1..N).each{|i|
  IO.foreach("gr.spk.#{i}"){|l|
    t, n = l.chomp.split
    indices.each_with_index{|m, j|
      files[j].puts "#{t.to_i*DT} #{i}" if n.to_i == m
    }
  }
}
        
indices.size.times{|i|
  files[i].close
}

indices.each{|i|
  ary = Array.new(T/Binsize)
  IO.foreach("gr.#{i}"){|l|
    t, n = l.chomp.split
    ary[(t.to_f/DT).to_i/Binsize] = ary[(t.to_f/DT).to_i/Binsize].to_i + 1
  }
  open("h.#{i}", "w"){|o|
    ary.each_with_index{|n, i|
      o.puts "#{DT*i.to_f*Binsize} #{n.to_f/(N.to_f)/(DT*Binsize)}"
    }
  }
}
