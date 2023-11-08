#include "bhns_header.h"

int BH_NS_Binary_Initial_Data(void *vp);
void bhns_export_id_generic(void *vp);

int BH_NS_Binary_Initial_Data(void *vp)
{
  /* if this is a generic ID reader call */
  if (strcmp_i(PgetsEZ("IDR_BHNS_export_id"),"generic"))
    bhns_export_id_generic(vp);
  else
    Error1(NO_OPTION);

  return EXIT_SUCCESS;
}

