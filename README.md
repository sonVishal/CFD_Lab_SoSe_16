# CFD_Lab_SoSe_16
Codes for the CFD lab at TUM for SoSe16 by Tord, Daniel and Vishal.

## Conventions

### Naming

|                       |             |
| ----------------------|-------------|
| **Macros, Constants**            |  MAX_BUFFER_SIZE  |
| **Struct and typedef's**    |  GtkWidget (-> camel case with first uppercase letter) |
| **Enum**              |  ETitleCase |
| **Enum Members**      |  ALL_CAPS |
| **Functions, default (operating on structs)**  |  my_function() (-> classic C style, underscores) |
| **private Functions (that are there, but shouldn't be called directly, or have obscure uses, e.g. 'helper functions')** |  _destroy_cache() (-> one or more underscores at the beginning) |
| **Trivial variables** |  i,x,n,f etc... |
| **Local variables**   |  myVariable (-> camel case with first lowercase letter) |
| **Global variables**  |  G_GLOBAL_VARIABLE |
