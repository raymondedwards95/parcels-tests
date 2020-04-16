""" Fast test for parcels """
from parcels import FieldSet, ParticleSet, AdvectionRK4, JITParticle
from datetime import timedelta

filenames = {"U": "GlobCurrent/*.nc",
             "V": "GlobCurrent/*.nc"}
var = {"U": "eastward_eulerian_current_velocity",
       "V": "northward_eulerian_current_velocity"}
dim = {"lat": "lat",
       "lon": "lon",
       "time": "time"}
fset = FieldSet.from_netcdf(filenames, variables=var, dimensions=dim, allow_time_extrapolation=True)

pset = ParticleSet(fieldset=fset, pclass=JITParticle, lon=[0], lat=[0])
print pset

pset.execute(
    AdvectionRK4,
    runtime=timedelta(days=150),
    dt=timedelta(minutes=10)
)

print pset
