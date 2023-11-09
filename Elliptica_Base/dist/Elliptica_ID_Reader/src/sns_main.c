#include "sns_header.h"

int Single_NS_Initial_Data(void *vp);
void sns_export_id_generic(void *vp);

int Single_NS_Initial_Data(void *vp)
{
  /* if this is a generic ID reader call */
  if (strcmp_i(PgetsEZ("IDR_SNS_export_id"),"generic"))
    sns_export_id_generic(vp);
  else
    Error1(NO_OPTION);

  return EXIT_SUCCESS;
}

