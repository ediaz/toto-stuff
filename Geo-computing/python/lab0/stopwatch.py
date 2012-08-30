#!/usr/bin/env python
import time as t


running = False 
time_ = 0.0
start_ = 0.0

def start():
    global running
    global time_ 
    global start_

    if not running:
        running = True
        start_ = t.time()

def stop(): 
    global running
    global time_ 
    global start_

    if (running) :
        time_ += t.time() - start_ 
        running = False


def reset(): 
    global running
    global time_ 
    global start_

    stop()
    time_ = 0.0

def restart():
    global running
    global time_ 
    global start_

    reset()
    start()

def time():
    global running
    global time_ 
    global start_

    if(running):
        return time_ + t.time() -start_
    else:
        return time_
