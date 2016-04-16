# CFD_Lab_SoSe_16
Codes for the CFD lab at TUM for SoSe16 by Tord, Daniel and Vishal.

## Conventions

### Naming

|                       |             |
| ----------------------|-------------|
| **Macros, Constants**            |  MAX_BUFFER_SIZE  |
| **Struct and Typedef's**    |  GtkWidget (-> camel case with first uppercase letter) |
| **Enum**              |  ETitleCase (-> prefix letter "E" indicates enum)  |
| **Enum Members**      |  ALL_CAPS |
| **Functions, default (operating on structs)**  |  my_function() (-> classic C style, underscores) |
| **Private Functions (that are there, but shouldn't be called directly, or have obscure uses, e.g. 'helper functions')** |  p_destroy_cache() (-> prefix "p_" to indicate private function (1)) |
| **Trivial variables** |  i,x,n,f etc... |
| **Local variables**   |  myVariable (-> camel case with first lowercase letter) |
| **Global variables**  |  G_GLOBAL_VARIABLE (-> use global variables with care (2)! prefix "G_" indicates global variable) |


(1) There are reserved C functions beginning with underscore only, therefore, "p_" is safer (no shadowing). http://c-faq.com/decl/namespace.html

(2) https://stackoverflow.com/questions/176118/when-is-it-ok-to-use-a-global-variable-in-c
 

### Line length (?)

Keep the line length to under 80 characters if possible. Then a split screen on
a laptop can show everything. 
