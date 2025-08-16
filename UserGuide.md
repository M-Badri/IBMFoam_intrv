# Documentation for IBMFoamDict

## 1) Global settings

**body_names** - *required* > List of body names (***body_name***)

**surface_threshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction)

**stepDEM** - *required* > Fraction of timestep used as a sub timestep for DEM calculation

**geometricD** - *optional* > Define empty direction for pseudo 2D simulation  
> *default value is taken from fvMesh i.e. empty patch is recognized.*

**time_scaling_factor** - *optional* > Factor for scaling simulation time without changing physical parameters

**save_simulation** - *required* > Should be particles recorded into bodiesInfo directory  
> Possible values: {*true*, *false*}

---

## 2) Interpolation schemes (near-body velocity reconstruction)

**interpolationSchemes** - *optional* > Interpolation settings used for velocity reconstruction near bodies

> **U** - *required* > Scheme for velocity interpolation inside a cell. Values according to Foam::interpolation
>> Possible values: {*cell*, *cellPoint*, *cellPointFace*}

> **method** - *required* > Method for interpolation point searching
>> Possible values: {*line*}

**gradients** - *optional* > Gradient interpolation control for auxiliary fields

---

## 3) Output controls

**out_setting** - *optional* > Output settings. When not provided, everything is outputted

> **basic** - *required* > Basic output
>> Possible values: {*true*, *false*}

> **iB** - *required* > Detail output for bodies
>> Possible values: {*true*, *false*}

> **DEM** - *required* > Output for DEM
>> Possible values: {*true*, *false*}

> **add_model** - *required* > Output for body addition
>> Possible values: {*true*, *false*}

**log_frequency** - *optional* > Interval (in timesteps) for writing status logs to file

---

## 4) DEM configuration

**DEM** - *required* > Input for DEM

> **materials** - *required* > Input for materials  
>> ***mat_name*** - *required* > Custom name of material Multiple blocks possible (see example)
>>> **Y** - *required* > Young's modulus  
>>> **nu** - *required* > Poisson ratio  
>>> **eps** - *required* > Viscoelastic damping constant  
>>> **mu** - *required* > Tangential force truncation  
>>> **adh_coef** - *required* > Adhesive force coefficient

> **face_adh** - *optional* > Truncation for adhesive for between materials  
>> ***name*** - *required* > Custom id name
>>> **materials** - *required* > Two names of affected materials  
>>> **value** - *required* > Value is used for both materials

> **coll_wall** - *required* > List of patches with which the bodies collide  
>> ***patchName*** - *required* > Name of the patch
>>> Possible values: ***mat_name***
>>> nVec: ***normal vector***
>>> pln_point: *** point located on plane***

**collision_damping_factor** - *optional* > Global damping factor applied to all collisions in DEM

---

## 5) Body block (repeat per body)

***body_name*** - *required* > Body name which correspond to value in **body_names**

> ***bodyType*** - *required* > Type of body
>> Possible values:
>> - *staticBody* > Body is static
>> - *prescribedTransBody* > Body has prescribed translational movement
>>> **velocity** - *required* > Translational velocity. (*in subdictionary*)
>> - *prescribedRotBody* > Body has prescribed rotational movement
>>> **axis** - *required* > Axis of rotation. (*in subdictionary*)  
>>> **omega** - *required* > Angular velocity . (*in subdictionary*)
>> - *prescribedTransRotBody* > Body has prescribed translational and rotational movement
>>> **velocity** - *required* > Translational velocity. (*in subdictionary*)  
>>> **axis** - *required* > Axis of rotation. (*in subdictionary*)  
>>> **omega** - *required* > Angular velocity . (*in subdictionary*)
>> - *prescribedTransFixedAxisRotBody* > Body has prescribed translational movement. Rotational movement has fixed axis.
>>> **velocity** - *required* > Translational velocity. (*in subdictionary*)  
>>> **axis** - *required* > Axis of rotation. (*in subdictionary*)
>> - *fullyCoupledBody* > Body is fully coupled with fluid phase
>>> **velocity** - *optional* > Translational velocity. (*in subdictionary*)

> **rho** - *required* > Density of the body

**thermal_expansion** - *optional* > Coefficient of thermal expansion applied to body geometry

> **U** - *required* > Dictionary for boundary condition
>> **BC** - *required* > Name of boundary condition
>>> Possible values: {*noSlip*}

> **material** - *required* > Material of the body
>> Possible values: ***mat_name***

> **b_geom** - *required* > Body geometry type
>> Possible values:
>> - *convex* > Body is defined according to STL file: ./constant/triSurface/***body_name***.stl. Solver is optimized for this type of geometry.
>> - *nonConvex* > Body is defined according to STL file: ./constant/triSurface/***body_name***.stl
>> - *sphere* > Body is defined by center and radius.
>> **sphere** - *required*  
>>> **startPosition** - *required* > Starting position of the center. (*in subdictionary*. May not be used in some **add_model**)  
>>> **radius** - *required* > Radius of the sphere. (*in subdictionary*)

> **update_torq** - *optional* > Should be rotational movement updated according to acting forces?
>> Possible values: {*true*, *false*}

> **interface_span** - *optional* > For sdBased algorithm. Best results with 1.0

> **startSynced** - *optional* > Should be body movement sync with fluid movement when added? Default false
>> Possible values: {*true*, *false*}

**color_id** - *optional* > Arbitrary color index for post-processing visualization

> **t_to_set_static** - *optional* > Number of timesteps to set body as static. Value -1 means never. Default -1

---

### Notes

- Define each body listed in `**body_names**` with its own ***body_name*** block under section 5.  
- For STL-based bodies (`*convex*` / `*nonConvex*`), ensure the file path and naming convention (`./constant/triSurface/***body_name***.stl`) are consistent with the declared `***body_name***`.  
- DEM sections can include multiple `***mat_name***` blocks to represent different materials; pairwise adhesive truncations can be specified via `**face_adh**`.  
- Optional keys like `thermal_expansion` and `color_id` are ignored if not specified.
