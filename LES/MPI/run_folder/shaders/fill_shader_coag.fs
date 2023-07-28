#version 330 core

out vec4 FragColor;

in vec3 ourColor;

vec3 output;

float val_sm = 1e-5;
float val_md = 1e-6;
float val_la = 1e-11;

float color;

void main()
{
    if (ourColor[2] > val_la)
    {
        color = ourColor[2]/(3.*val_la);
        if (color > 1.) color = 1.;
        output = vec3(color, 0.0, 1. - color);
    }
    else if (ourColor[1] > val_md)
    {
        color = ourColor[1]/(3.*val_md);
        if (color > 1.) color = 1.;
        output = vec3(color, color, 1. - color);
    }
    else if (ourColor[0] > val_sm)
    {
        color = ourColor[0]/(3.*val_sm);
        if (color > 1.) color = 1.;
        output = vec3(0.0, color, 1. - color);
    }    
    else
    {
        output = vec3(1.,1.,1.);
    }

    FragColor = vec4(output, 1.0);
    
}
