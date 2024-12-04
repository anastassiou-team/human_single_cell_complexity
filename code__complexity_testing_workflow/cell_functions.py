"""Additional axon processing options for all-active models"""


from neuron import h
from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor
from bmtk.simulator.bionet.default_setters.cell_models import set_params_allactive


def fix_axon_noaxon(hobj):

    for sec in hobj.axon:
        h.delete_section(sec=sec)
    h.define_shape()


def csmc_allactive_noaxon(hobj, cell, dynamics_params):
    fix_axon_noaxon(hobj) 
    set_params_allactive(hobj, dynamics_params)
    return hobj

add_cell_processor(csmc_allactive_noaxon)
