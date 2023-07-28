#version 330 core

out vec4 FragColor;

in vec3 ourColor;

vec3 output;

float val_sm = 1e-9;
vec3 color_sm = vec3(82./255., 219./255., 1.);
float val_md = 1e-17;
vec3 color_md = vec3(1., 173./255., 47./255.);
float val_la = 1e-24;
vec3 color_la = vec3(127./255., 0., 1.);

vec3 color;
float val;

float main_color;

void main()
{
    if (ourColor[0] > val_sm)
    {
        main_color = ourColor[0];
        val = val_sm*10;
        color = color_sm;
    }    
    if (ourColor[1] > val_md)
    {
        main_color = ourColor[1];
        val = val_md*10;
        color = color_md;
    }
    if (ourColor[2] > val_la)
    {
        main_color = ourColor[2];
        val = val_la*10;
        color = color_la;
    }

    if (ourColor[0] <= val_sm && ourColor[1] <= val_md && ourColor[2] <= val_la)
    {
        main_color = ourColor[0];
        val = -1;
        color = vec3(1.,1.,1.);
    }


    //if (main_color < val)
    //{
    //    output = vec3(0.0, (1./val)*main_color, (1./val)*(val-main_color));
    //}
    //else
    //{
        output = color;
    //}

    FragColor = vec4(output, 1.0);
    
}
