#!/usr/bin/env python

import time as t 
class stopwatch:

    def __init__(self):
        self.running= False
        self.time_  = 0.0 
        self.start_ = 0.0

    def start(self):
        global t

        if not self.running:
            self.running = True
            self.start_ = t.time()
    
    def stop(self): 
        global t
    
        if (self.running) :
            self.time_ += t.time() - self.start_ 
            self.running = False
    
    def reset(self): 
        global t
        self.stop()
        self.time_ = 0.0
    
    def restart(self):
        global t
        self.reset()
        self.start()
    
    def time(self):
        global t
        if(self.running):
            return self.time_ + t.time() -self.start_
        else:
            return self.time_
