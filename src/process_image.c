#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{   
    // implement clamp
    if (x > im.w) x = im.w-1;
    if (x < 0)    x = 0;

    if (y > im.h) y = im.h-1;
    if (y < 0)    y = 0;

    if (c > im.c) c = im.c-1;
    if (c < 0)    c = 0;
    // Channel, row, col ordering
    // channel_of_interest * width * height = the start of the channel we want 
    // row_of_interest * width = the start of the row we want
    int index = c*im.w*im.h + y*im.w + x;
    return im.data[index];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if ((x>im.w) || (x<0) || (y>im.h) || (y<0) || (c>im.c) || (c<0)){
        return;
    }
    
    int index = c*im.w*im.h + y*im.w + x;
    im.data[index] = v;

}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data,im.data,__SIZEOF_FLOAT__*im.w*im.h*im.c);
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);

    int ii,jj;
    for(jj=0; jj<im.h; ++jj){
        for(ii=0;ii<im.w;++ii){
            int dst_idx = jj*im.w + ii;
            int src_idx = jj*im.w + ii;
            int size = im.w*im.h;
            gray.data[dst_idx] = 0.299*im.data[src_idx] + 
                                 0.587*im.data[src_idx+size] + 
                                 0.114*im.data[src_idx+2*size];
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    int ii,jj;
    for(jj=0; jj<im.h; ++jj){
        for(ii=0;ii<im.w;++ii){
            int idx = c*im.w*im.h + jj*im.w + ii;
            im.data[idx] += v;
        }
    }
}

void clamp_image(image im)
{
    for(int idx=0; idx<im.w*im.h*im.c; ++idx){
        if(im.data[idx] > 1.0) im.data[idx] = 1.0;
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    int img_size = im.w*im.h;

    float H,S,V,C,m;
    for(int jj=0; jj<im.h; ++jj){
        for(int ii=0; ii<im.w; ++ii){
            int r_idx =               jj*im.w + ii;
            int g_idx =    img_size + jj*im.w + ii;
            int b_idx =  2*img_size + jj*im.w + ii;

            V = three_way_max(im.data[r_idx],im.data[g_idx],im.data[b_idx]);
            m = three_way_min(im.data[r_idx],im.data[g_idx],im.data[b_idx]);
            C = V-m;
            
            //if(im.data[r_idx]==0.0 && im.data[g_idx]==0.0 && im.data[b_idx]==0.0){
            if(C==0.0 || V == 0.0){
                S = 0.0;
            }else{
                S = C/V;
            }

            if (C==0.0){
                H = 0.0;
            }else if (V==im.data[r_idx])
            {
                H = (im.data[g_idx]-im.data[b_idx])/C;
            }else if (V==im.data[g_idx])
            {
                H = (im.data[b_idx]-im.data[r_idx])/C+2;
            }else //(V==im.data[b_idx])
            {
                H = (im.data[r_idx]-im.data[g_idx])/C+4;
            }

            H = H < 0 ? H/6.0+1.0 : H/6.0;
            im.data[r_idx] = H;
            im.data[g_idx] = S;
            im.data[b_idx] = V;
        }
    }
}

void hsv_to_rgb(image im)
{
    ///https://www.rapidtables.com/convert/color/hsv-to-rgb.html
    int img_size = im.w*im.h;

    float H,S,V,C,m,X;
    //int H,X;
    float R,G,B;
    for(int jj=0; jj<im.h; ++jj){
        for(int ii=0; ii<im.w; ++ii){
            int h_idx =               jj*im.w + ii;
            int s_idx =    img_size + jj*im.w + ii;
            int v_idx =  2*img_size + jj*im.w + ii;

            H = im.data[h_idx] * 360;
            S = im.data[s_idx];
            V = im.data[v_idx];

            C = S*V;
            m = V - C;
            X = C*(1-fabs(fmod((H/60),2)-1));

            if(H<60){
                R=C;
                G=X;
                B=0;
            }else if (H<120)  
            {
                R=X;
                G=C;
                B=0;
            }else if (H<180)  
            {
                R=0;
                G=C;
                B=X;
            }else if (H<240)  
            {
                R=0;
                G=X;
                B=C;
            }else if (H<300)  
            {
                R=X;
                G=0;
                B=C;
            }else
            {
                R=C;
                G=0;
                B=X;
            }
            
            im.data[h_idx] = R+m;
            im.data[s_idx] = G+m;
            im.data[v_idx] = B+m;
        }
    }
}
