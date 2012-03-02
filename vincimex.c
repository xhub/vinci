


#include "mex.h"
#include "vinci.h"

#define VINCI_VERSION 1.0.5


/* functions to interface with vinci */

void vertMatlabToVinci(const mxArray *mat)
{
	T_Vertex *v;
	double * vArray;
	unsigned int i;
	unsigned int j;
	G_n = mxGetM(mat);
	G_d = mxGetN(mat);
	vArray = mxGetPr(mat);
	G_Vertices = create_empty_set();
	for (i = 0; i < G_n; i++)
	{
		v = create_vertex ();
		v->no = i;
		for (j = 0; j < G_d; j++)
		{
			v->coords[j] = vArray[j*G_n+i];
		}
		add_element(&G_Vertices, v);
	}
}

void HMatlabToVinci(const mxArray *matH, const mxArray *vecK)
{
	double *H;
	double *K;
	unsigned int i;
	unsigned int j;
	G_m = mxGetM(matH);
	G_d = mxGetN(matH);
	create_hyperplanes();
	H = mxGetPr(matH);
	K = mxGetPr(vecK);
	for (i = 0; i < G_m; i++)
	{
		G_Hyperplanes[i][G_d] = K[i];
			for (j = 0; j < G_d; j++)
			{
				G_Hyperplanes[i][j] = - H[j*G_m+i];
			}
	}
}
/* Function needed by matlab */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray *matVertices;
	const mxArray *matH;
	const mxArray *vecK;
	char *input_buf;
	boolean Vrep = 0;
	boolean Hrep = 0;
	boolean methodOK;
	rational *volume;
	int method, bufLen;
	char *inputBuf;

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	volume = mxGetPr(plhs[0]);
	if (nrhs == 4)
	{
		/* we are using a V representation */
		/* get the vertices */
		matVertices = prhs[1];
		vertMatlabToVinci(matVertices);
		matH = prhs[2];
		vecK = prhs[3];
		HMatlabToVinci(matH, vecK);
	}
	else
	{
		/* Wrong */
		mexErrMsgTxt("You have to give 4 arguments");
		return;
	}
	/* First arg must be the method */
	if (mxIsChar(prhs[0]) && (mxGetM(prhs[0]) == 1))
	{
		bufLen = mxGetN(prhs[0]) + 1;
		inputBuf = mxCalloc(bufLen, sizeof(char));
		mxGetString(prhs[0], inputBuf, bufLen);
		/* see which method we want to use */
		if (strcmp(inputBuf,"rch\0") == 0)
			method = RCH;
		else if (strcmp(inputBuf,"hot\0") == 0)
			method = HOT;
		else if (strcmp(inputBuf,"lawnd\0") == 0)
			method = LAWND;
		else if (strcmp(inputBuf,"lawd\0") == 0)
			method = LAWD;
		else if (strcmp(inputBuf,"rlass\0") == 0)
			method = RLASS;
		else if (strcmp(inputBuf,"lrs\0") == 0)
			method = LRS;
		else
			mexErrMsgTxt("To be implemented");

		switch (method)
			{
/*				case RCH:
					volChV(volume);
					break;*/
				case HOT:
					volOrthoV(volume);
					break;
/*				case LAWND:
					volLawrenceV(volume);
					break;
				case LAWD:
					volLawrenceLrsV(volume);
					break;
				case RLASS:
					volLasserreV(volume);
					break;
				case LRS:
					volLrsV(volume);
					break;*/
			}
/*		if (Hrep == 1)
		{
			switch (method)
			{
				case RCH:
					volChH(volume);
					break;
				case HOT:
					volOrthoH(volume);
					break;
				case LAWND:
					volLawrenceH(volume);
					break;
				case LAWD:
					volLawrenceLrsH(volume);
					break;
				case RLASS:
					volLasserreH(volume);
					break;
				case LRS:
					volLrsH(volume);
					break;
			}
		}*/
	}
	else
	{
		mexErrMsgTxt("The first argument must be a string : it is the name \
				we want to use.\n Leave it empty if you want to use the autodetect method");
	}
}
