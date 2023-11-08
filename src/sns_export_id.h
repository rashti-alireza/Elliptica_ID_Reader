#include "sns_header.h"
#include "id_reader_lib.h"
#include "elliptica_id_reader_lib.h"

/* prefix parameters came from evo codes, should be lower case  */
#define EVO_ "evo_"

#define STR_LEN_MAX (9999)

void sns_export_id_generic(void *vp);
void sns_set_evo_fields_generic(Grid_T *const grid);

