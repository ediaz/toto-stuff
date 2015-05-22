module Stopwatch
#
#
# Stopwatch module,
# to get the object, user must call
# Stopwatch.init() which returns an StopWatch object
# which can be passed later to stop, reset, restart
# and get_time() functions
#
#
type StopWatch
  running::Bool
  time::Float32
  start::Float32   
end

function init()
  return StopWatch(true,0.0,time())
end

function stop(sw:: StopWatch)
  if sw.running 
    sw.time = sw.time +time() - sw.start
    sw.running = false
  end
end

function reset(sw:: StopWatch)
  stop(sw)
  sw.time = 0
end

function restart(sw:: StopWatch)
  sw.running = true
  sw.time = 0.0
  sw.start = time() 
end


function get_time(sw:: StopWatch)
  if sw.running
    t = sw.time - sw.start + time()
  else
    t = sw.time
  end
  return t
end

end
