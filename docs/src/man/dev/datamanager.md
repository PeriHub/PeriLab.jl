# Datamanager
The datamanager is a central part of PeriLab. You can create fields as you need them. 
## Create constant node fields
Constant node fields are fields with the length of the number of nodes. It return one vector of type Type_of_variable. 

```julia
field = datamanager.create_constant_node_field(Fieldname::String, Type_of_variable::Type, Degree_of_freedom::Int64)
```

You can get constant fields  
```julia
field = datamanager.get_field(Fieldname::String)
```

## Create non-constant node fields
Non-constant node fields are fields with the length of the number of nodes. Each node could have any degree of freedom. 
You will get a field with N and N+1 created it by this command

```julia
fieldN, fieldNP1 = datamanager.create_node_field(Fieldname::String, Type_of_variable::Type, Degree_of_freedom::Int64)
```
!!! info "Switch values" The non-constant fields are switched automatically at the end of each time integration step. At the beginning of the next step NP1 is zero, wheras N is the NP1 from the previous step.

You can get non-constant fields as   
```julia
fieldN = datamanager.get_field(Fieldname::String, "N")
fieldNP1 = datamanager.get_field(Fieldname::String, "NP1")
fieldN = datamanager.get_field(Fieldname*"N"::String)
fieldNP1 = datamanager.get_field(Fieldname*"NP1"::String)
```
## Create constant bond fields
Constant bond fields are fields with the length of the number of nodes. Each node has a vector of lenght number of neighbors or bonds with a defined degree of freedom. It return one vector of type Type_of_variable.

```julia
field = datamanager.create_constant_node_field(Fieldname::String, Type_of_variable::Type, Degree_of_freedom::Int64)
```

You can get constant fields  
```julia
field = datamanager.get_field(Fieldname::String)
```

## Create non-constant bond fields
Non-constant node fields are fields with the length of the number of nodes. Each node has a vector of lenght number of neighbors or bonds with a defined degree of freedom. 
You will get a field of type Type_of_variable with N and N+1 created it by this command

```julia
fieldN, fieldNP1 = datamanager.create_bond_field(Fieldname::String, Type_of_variable::Type, Degree_of_freedom::Int64)
```
!!! info "Switch values" The non-constant fields are switched automatically at the end of each time integration step. At the beginning of the next step NP1 is zero, wheras N is the NP1 from the previous step.

You can get non-constant fields as   
```julia
fieldN = datamanager.get_field(Fieldname::String, "N")
fieldNP1 = datamanager.get_field(Fieldname::String, "NP1")
fieldN = datamanager.get_field(Fieldname*"N"::String)
fieldNP1 = datamanager.get_field(Fieldname*"NP1"::String)
```

## Other Options
For node and bond fields (constant and non-constant) the following options are possible.

**Matrix style**

You can switch between vector and matrix. If you give the keyword Matrix_or_Vector = "Vector" you get vector of length degree of freedom for each node or bond. If Matrix_or_Vector = "Matrix" each node or bond gets matrix dof $\times$ dof
```julia
datamanager.create_node_field(Fieldname::String, Type_of_variable::Type, Matrix_or_Vector::String, Degree_of_freedom::Int64)
```
**Pre-defined Values**

You can add an optional value at the end. The default is zero.

```julia
datamanager.create_node_field(Fieldname::String, Type_of_variable::Type, Degree_of_freedom::Int64, Value::Type_of_variable)
```

## Create free size fields
Free size fields can be constant or non-constant. The size is defined by the user and must be a tuple. The field dimension is chosen by the user.
```julia
 field = datamanager.create_constant_free_size_field(Fieldname::String, Type_of_variable::Type, size::Tuple)
 fieldN, fieldNP1 = datamanager.create_free_size_field(Fieldname::String, Type_of_variable::Type, size::Tuple)
```
!!! info "free size example" 
    A constant node field with a matrix $3\times3$ per node can be defined in a free size field by


```julia
field = datamanager.create_constant_free_size_field("Example", Float64, (number_of_nodes, 3, 3))
```