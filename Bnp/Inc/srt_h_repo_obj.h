
#ifndef SRT_H_REPO_OBJ_H
#define SRT_H_REPO_OBJ_H

typedef SrtListHdr SrtRepoObj;

/* Build a SrtRepoObj from a vector of repos */
Err srt_f_repo_obj_init(
			Date		today,
			double		**repos,
			int			num_cols,
			int			num_rows,
			char		*repo_curve_name,
			SrtRepoObj	**repo_obj);

Err srt_f_repo_obj_repo(SrtRepoObj *repo_obj,
						double start_date,
						double *repo);

Err srt_f_repo_obj_create(
			Date		today,
			char		*repo_crv_name,
			SrtRepoObj	**repo_obj);

Err srt_f_repo_obj_delete(SrtRepoObj	**repo_obj);

Err srt_f_repo_obj_insert_point(SrtRepoObj	*repo_obj,Ddate date,double repo);

Err srt_f_repo_obj_extracttoday(SrtRepoObj *repo_obj, Ddate *today);


#endif
