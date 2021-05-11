import numpy as np, matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import ICRS
from astropy.time import Time
import sys
import yaml
from astroplan import Observer
import astroplan
from astroplan import Constraint, is_observable, min_best_rescale, max_best_rescale
from astroplan.constraints import _get_altaz
from astroplan import FixedTarget
from astroplan import ObservingBlock
from astroplan.constraints import TimeConstraint, AirmassConstraint
from astroplan.scheduling import Transitioner
from astroplan.scheduling import SequentialScheduler, PriorityScheduler
from astroplan.scheduling import Schedule

class AzimuthConstraint(Constraint):
    """
    Constrain the azimuth 
    """
    def __init__(self, min=None, max=None, boolean_constraint=True):
        """
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable azimuth. `None`
            indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable azimuth. 
        """
        self.min = min if min is not None else 0*u.deg
        self.max = max if max is not None else 90*u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times, observer, targets):

        cached_altaz = _get_altaz(times, observer, targets)
        az = cached_altaz['altaz'].az
        fixed_az = np.arccos(np.abs(np.cos(az)))
        
        if self.boolean_constraint:
            lowermask = self.min <= fixed_az
            uppermask = fixed_az <= self.max
            return lowermask & uppermask
        else:
            return min_best_rescale(fixed_az.value*180./np.pi, self.min.value, self.max.value)
        
from astropy.table import Table

def to_table(sch, show_transitions=True, show_unused=False):
    # TODO: allow different coordinate types                                                                  
    target_names = []
    start_times = []
    end_times = []
    durations = []
    ra = []
    dec = []
    config = []
    for slot in sch.slots:
        if hasattr(slot.block, 'target'):
            start_times.append(slot.start.iso)
            end_times.append(slot.end.iso)
            durations.append(slot.duration.to(u.minute).value)
            target_names.append(slot.block.target.name)
            ra.append(slot.block.target.ra)
            dec.append(slot.block.target.dec)
            config.append(slot.block.configuration)
        elif show_transitions and slot.block:
            start_times.append(slot.start.iso)
            end_times.append(slot.end.iso)
            durations.append(slot.duration.to(u.minute).value)
            target_names.append('TransitionBlock')
            ra.append('')
            dec.append('')
            changes = list(slot.block.components.keys())
            if 'slew_time' in changes:
                changes.remove('slew_time')
            config.append(changes)
        
    return Table([target_names, start_times, end_times, durations,(ra), (dec), config],
                 names=('target', 'start time (UTC)', 'end time (UTC)',
                        'duration (minutes)', 'ra', 'dec', 'configuration'))

def sortFunc(e):
    return e['time'].timestamp()

# to convert to series of actions...
# use transitions to deal with moves
# between transitions, identify start/stop times
# 
def define_actions(tab,sched_start,recording=False):
    
    transition_time = 290*u.second
    transition_time2 = 10*u.second
    start_time = 60*u.second
    stop_time = 180*u.second
    nominal_obs_time = 3600*u.second
    
    actions = []
    
    for ct in np.arange(len(tab)):
        
        if ct==0:
            delta_obs = Time(tab[ct]['end time (UTC)']) - Time(sched_start) - transition_time
            tmove = Time(sched_start)
        else:
            delta_obs = Time(tab[ct]['end time (UTC)']) - Time(tab[ct-1]['end time (UTC)']) - transition_time
            tmove = Time(tab[ct-1]['end time (UTC)'])
            
        psr=False
        if tab[ct]['configuration']['dm'] is not None:
            
            tm = Time(tab[ct]['start time (UTC)']) + 750*u.second
            if recording is True:
                actions.append({'time':tm.datetime, 'cmd':'record','val':'300-'+tab[ct]['target']+'-'})
            
        n_block = (delta_obs.to_value('s')/(start_time + stop_time + nominal_obs_time)).value
        full_blocks = int(np.floor(n_block))
        part_block = n_block - np.floor(n_block)
        
        # move action
        actions.append({'time':tmove.datetime, 'cmd':'move','val':90.-(37.23-tab[ct]['dec'].deg),'config':tab[ct]['configuration']})
        tmove += transition_time
        actions.append({'time':tmove.datetime, 'cmd':'move','val':90.-(37.23-tab[ct]['dec'].deg),'config':tab[ct]['configuration']})
        tmove += transition_time2
        
        # psr?
        psr=False
        if tab[ct]['configuration']['dm'] is not None:
            psr=True
        
        # start/stop actions
        for i in np.arange(full_blocks):
            
            actions.append({'time':tmove.datetime, 'cmd':'start','val':0})
            tmove += start_time + nominal_obs_time
            
            if i<full_blocks-1:
                actions.append({'time':tmove.datetime, 'cmd':'stop','val':0})
                tmove += stop_time
                
            elif i==full_blocks-1:
                if part_block<=0.5:
                    extra_obs_time = part_block*(start_time + stop_time + nominal_obs_time)
                    tmove += extra_obs_time
                    
                    actions.append({'time':tmove.datetime, 'cmd':'stop','val':0})
                    tmove += stop_time
                    
                else:
                    actions.append({'time':tmove.datetime, 'cmd':'stop','val':0})
                    tmove += stop_time
                    
                    
            
        if part_block>0.5:
            
            obs_time = part_block*(start_time + stop_time + nominal_obs_time) - start_time - stop_time
            
            actions.append({'time':tmove.datetime, 'cmd':'start','val':0})
            tmove += start_time + obs_time
            actions.append({'time':tmove.datetime, 'cmd':'stop','val':0})
            tmove += stop_time
            
    a = sorted(actions, key = lambda i: i['time'].timestamp())
    return a


# read in sources
def read_srcs(fl):
    with open(fl, 'r') as stream:
        try:
            srcs = yaml.load(stream)
            return(srcs)
        except yaml.YAMLError as exc:
            print('cannot open yaml file')