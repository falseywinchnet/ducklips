// Preserved user configuration and its below-mouth extension.
export const userReference = {
  "version": 2,
  "pressure": 570,
  "wander": 0.06,
  "monitor": 0.25,
  "blocks": [],
  "bands": [
    {
      "name": "01 \u00b7 Return chamber",
      "length": 0.05,
      "radius": 0.009,
      "material": "soft",
      "mass": 1.05,
      "elasticity": 1.85,
      "offset": -0.0014740847632998513,
      "anchor": 0.5401597282112418,
      "support": 0.85
    },
    {
      "name": "02 \u00b7 Vocal band",
      "length": 0.016,
      "radius": 0.00032,
      "material": "fold",
      "mass": 1.85,
      "elasticity": 3.4,
      "offset": 0.00027209617978050596,
      "anchor": 0.549436346277014,
      "support": 0.2
    },
    {
      "name": "03 \u00b7 Lower tract",
      "length": 0.035,
      "radius": 0.007,
      "material": "soft",
      "mass": 3.5,
      "elasticity": 1.5,
      "offset": -0.004313555762881324,
      "anchor": 0.5029170886660774,
      "support": 0.2
    },
    {
      "name": "04 \u00b7 Oral body",
      "length": 0.045,
      "radius": 0.011,
      "material": "firm",
      "mass": 4.7,
      "elasticity": 0.2,
      "offset": -0.008342391967773437,
      "anchor": 0.5373209353580447,
      "support": 3.2
    },
    {
      "name": "05 \u00b7 Front cavity",
      "length": 0.035,
      "radius": 0.008,
      "material": "fold",
      "mass": 0.2,
      "elasticity": 3,
      "offset": -0.0063339342389787945,
      "anchor": 0.5297026512131847,
      "support": 0.2
    },
    {
      "name": "06 \u00b7 Outlet",
      "length": 0.025,
      "radius": 0.006,
      "material": "soft",
      "mass": 1,
      "elasticity": 2.1,
      "offset": 0.002478314354306176,
      "anchor": 0.5457723404525157,
      "support": 1
    }
  ]
};
export function referenceDesign(){return JSON.parse(JSON.stringify(userReference));}
export function belowMouthDesign(){
 const d=referenceDesign(),[inlet,cord,soft,firm,adapter,outlet]=d.bands;
 const radius=soft.radius+soft.offset;
 const ring=(b,name,role,changes={})=>({...b,name,role,points:8,...changes});
 const narrow=b=>({offset:radius-b.radius});
 const lowAdapter={...narrow(adapter),elasticity:3};
 return {...d,version:3,system:'below-mouth',pressure:850,elasticityScale:1,elasticityVariation:.10,noiseSeed:73021,blocks:[],bands:[
  ring(inlet,'01 · Inlet expansion','inlet',{offset:2*radius-inlet.radius,mass:1,elasticity:1,support:1}),
  ring(adapter,'02 · Lower elastic adapter','adapter',{...lowAdapter,length:adapter.length/2,points:4}),
  ring(firm,'03 · Lower firm support','firm',{...narrow(firm),length:firm.length/2,points:4}),
  ring(soft,'04 · Lower heavy tissue','soft',{length:soft.length/2,points:4}),
  ring(cord,'05 · Vocal cord','cord',{elasticity:1.5}),
  ring(soft,'06 · Upper heavy tissue','soft'),
  ring(firm,'07 · Upper firm support','firm',narrow(firm)),
  ring(adapter,'08 · Upper elastic adapter','adapter',lowAdapter),
  ring(outlet,'09 · Throat expansion','outlet',{offset:2*radius-outlet.radius,mass:1,elasticity:1,support:1})
 ]};
}
