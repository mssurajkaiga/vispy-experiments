# -*- coding: utf-8 -*-
# Copyright (c) 2014, Vispy Development Team.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

from __future__ import division

from ._default_app import default_app
from ..util.event import Event, EmitterGroup
from ..util.ptime import time as precision_time
from .base import BaseTimerBackend as TimerBackend  # noqa


class Timer(object):

    """Timer used to schedule events in the future or on a repeating schedule.
    """

    def __init__(
            self,
            interval=0.0,
            connect=None,
            iterations=-1,
            start=False,
            app=None):
        self.events = EmitterGroup(source=self,
                                   start=Event,
                                   stop=Event,
                                   timeout=Event)
        #self.connect = self.events.timeout.connect
        #self.disconnect = self.events.timeout.disconnect

        # Get app instance and make sure that it has an associated backend
        self._app = default_app if app is None else app
        self._app.use()

        # Instantiate the backed with the right class
        self._backend = self._app.backend_module.TimerBackend(self)

        self._interval = interval
        self._running = False
        self._last_emit_time = None
        self.iter_count = 0
        self.max_iterations = iterations
        if connect is not None:
            self.connect(connect)
        if start:
            self.start()

    @property
    def app(self):
        """ The vispy Application instance on which this Timer is based.
        """
        return self._app

    @property
    def interval(self):
        return self._interval

    @interval.setter
    def interval(self, val):
        self._interval = val
        if self.running:
            self.stop()
            self.start()

    @property
    def running(self):
        return self._running

    def start(self, interval=None, iterations=None):
        """Start the timer.

        A timeout event will be generated every *interval* seconds.
        If *interval* is None, then self.interval will be used.

        If *iterations* is specified, the timer will stop after
        emitting that number of events. If unspecified, then
        the previous value of self.iterations will be used. If the value is
        negative, then the timer will continue running until stop() is called.
        """
        self.iter_count = 0
        if interval is not None:
            self.interval = interval
        if iterations is not None:
            self.max_iterations = iterations
        self._backend._vispy_start(self.interval)
        self._running = True
        self._last_emit_time = None
        self.events.start(type='timer_start')

    def stop(self):
        """Stop the timer."""
        self._backend._vispy_stop()
        self._running = False
        self.events.stop(type='timer_stop')

    # use timer.app.run() and .quit() instead.
    # def run_event_loop(self):
        #"""Execute the event loop for this Timer's backend.
        #"""
        # return self._backend._vispy_run()

    # def quit_event_loop(self):
        #"""Exit the event loop for this Timer's backend.
        #"""
        # return self._backend._vispy_quit()

    @property
    def native(self):
        """ The native timer on which this Timer is based.
        """
        return self._backend._vispy_get_native_timer()

    def _timeout(self, *args):
        # called when the backend timer has triggered.
        if not self.running:
            return
        if self.max_iterations >= 0 and self.iter_count >= self.max_iterations:
            self.stop()
            return

        # compute dt since last event
        now = precision_time()
        if self._last_emit_time is None:
            dt = None
        else:
            dt = now - self._last_emit_time
        self._last_emit_time = now

        self.events.timeout(
            type='timer_timeout',
            iteration=self.iter_count,
            dt=dt)
        self.iter_count += 1

    def connect(self, callback):
        """ Alias for self.events.timeout.connect() """
        return self.events.timeout.connect(callback)

    def disconnect(self, callback=None):
        """ Alias for self.events.timeout.disconnect() """
        return self.events.timeout.disconnect(callback)
