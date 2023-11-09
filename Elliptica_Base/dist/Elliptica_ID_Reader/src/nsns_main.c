#include "nsns_header.h"

int NS_NS_Binary_Initial_Data(void *vp);
void nsns_export_id_generic(void *vp);

int NS_NS_Binary_Initial_Data(void *vp)
{
  /* if this is a generic ID reader call */
  if (strcmp_i(PgetsEZ("IDR_NSNS_export_id"),"generic"))
    nsns_export_id_generic(vp);
  else
    Error1(NO_OPTION);

  return EXIT_SUCCESS;
}

