#!/usr/bin/ruby

def psth(trial, cell_number)
  sum = Array.new(2000,0)
  ((trial-9)..trial).each{|k|
    name = "pkjvnio.spk.#{k}"
    open(name){|f|
      while line = f.gets do
        t, n = line.split
        t = t.to_i
        n = n.to_i
        if n == cell_number
          sum[t] = sum[t].to_i + 1
        end
      end
    }
  }
  sum2 = Array.new(20,0)
  (0...20).each{|i|
    n = 0
    (0...100).each{|j|
      n = n + sum[j+100*i]/1000.0
    }
    sum2[i] = n
  }
  open("#{trial}_#{cell_number}.dat", "w"){|o|
    sum2.each_with_index{|n, i|
      o.puts "#{i*0.1} #{n*1000}"
    }
  }
end


def main
  trial = ARGV[0].to_i
  cell_number = if ARGV.size == 2 then ARGV[1].to_i else 16 end
  psth(trial, cell_number)
end

main if __FILE__ == $0
