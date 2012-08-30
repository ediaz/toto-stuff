package lab0;

/**
 * A timer that works like a common stopwatch.
 * @author Dave Hale, Colorado School of Mines
 * @version 2004.11.02
 */
public class Stopwatch {

  public void start() {
    if (!_running) {
      _running = true;
      _start = System.nanoTime();
    }
  }

  public void stop() {
    if (_running) {
      _time += System.nanoTime()-_start;
      _running = false;
    }
  }

  public void reset() {
    stop();
    _time = 0;
  }

  public void restart() {
    reset();
    start();
  }

  public double time() {
    if (_running) {
      return 1.0e-9*(_time+(System.nanoTime()-_start));
    } else {
      return 1.0e-9*_time;
    }
  }

  private boolean _running;
  private long _start;
  private long _time;
}
