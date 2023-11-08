#include "bhbh_header.h"

int BH_BH_Binary_Initial_Data(void *vp);
void bhbh_export_id_generic(void *vp);

int BH_BH_Binary_Initial_Data(void *vp)
{
  /* if this is a generic ID reader call */
  if (strcmp_i(PgetsEZ("IDR_BHBH_export_id"),"generic"))
    bhbh_export_id_generic(vp);
  else
    Error1(NO_OPTION);

  return EXIT_SUCCESS;
}

