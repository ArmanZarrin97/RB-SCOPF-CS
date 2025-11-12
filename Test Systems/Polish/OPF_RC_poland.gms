$Title optimal power flow model with voltages expressed in rectangular coordinates
$OnText

Authors : Florin Capitanescu (florin.capitanescu@uni.lu) 
Date : XX October 2013
(c) 2013 University of Luxembourg,
         Interdisciplinary Centre for Security Reliability and Trust (SnT)
All rights reserved
GPL version X.XX to be found at: www.xxx.com

Comments: although the model is nonlinear due to products status of breakers times power flows, 
          the model has good relaxation, and converges acceptably fast 

$OffText

$offlisting
$offsymxref
$offsymlist

option

* Equations listed per block
         limrow = 0
* Variables listed per block
       , limcol = 0
* Solver's solution output printed
       , solprint = on
* Solver's system output printed
       , sysout = off
* Use CPLEX
*       , QCP = Cplex
       ;

Sets

         k       Contingencies / 1 /

         i       Nodes / 1*2746 /
         ichar   Nodes characteristics / "Vmin","Vnom","Vmax","BusN","e0","Vim0" /

         f       Branches / LINE0001*LINE3107 /
         bchar   Branches characteristcs / "From","To","g","b","bsh","Snom","BStatus0","dBStatus" /

         t       Transformers / T001*T172 /
         tchar   Transformers characteristcs / "From","To","g","b","bsh","ratio","rmin","rmax","Snom","BStatus0" /

         l       Loads /  L17*L2746 /
         lchar   Loads characteristics / "Bus","Pload","QLoad","ShFr","ShFrMin","ShFrMax","ShFrOK" /

         g       Generators / G17*G2740 /
         gchar   Generators characteristics / "Bus","P0","Q0","Vgen","Snom","Pmin","Pmax","Qmin","Qmax","genP","genQ","c0","c1","c2","BStatus" /

*         s       Shunts / S1*S14 /
*         schar   Shunts characteristics / "Bus","bsh0","bshmin","bshmax" /

         ;

         alias(j,i);

Scalars

         NL Number of lines                                  / 3107 /
         VR Voltage equation relaxation  (0 = No relaxation) / 0 /

         ;

$include C:\opf\spooles\linsol\drivers\code_gams\bus_poland.gms
$include C:\opf\spooles\linsol\drivers\code_gams\line_poland.gms
$include C:\opf\spooles\linsol\drivers\code_gams\load_poland.gms
$include C:\opf\spooles\linsol\drivers\code_gams\gen_poland.gms
*$include C:\opf\spooles\linsol\drivers\code_gams\shunt118.gms
$include C:\opf\spooles\linsol\drivers\code_gams\trfo_poland.gms

Variables

         Pg(g,k) Active power generation [pu]
         Qg(g,k) Reactive power generation [pu]
        
*         Pij(f,k) [pu]
*         Qij(f,k) [pu]

         ratio(t,k) [pu] transformer ratio
*         bshunt(s,k) [pu] shunt susceptance

         e(i,k) [pu] real part of bus voltage
         Vim(i,k) [pu] imaginary part of bus voltage
*         voltage(i,k) [pu] voltage magnitude

*         BStatus(f,k) Status of a branch [On-Off]
*         BStatuschange(f,k) Change of branch status [binary] 

         Cost Objective function value

         positive variables e,ratio ;
*         Binary variables BStatuschange ;

         Pg.l(g,k) = gendata(g,'P0');
         Qg.l(g,k) = gendata(g,'Q0');

         e.l(i,k) = idata(i,'e0');
         Vim.l(i,k) = idata(i,'Vim0');
         Vim.fx('28',k) = 0 ;

*         BStatuschange.l(f,k)=0;
*         BStatus.l(f,k) = bdata(f,'BStatus0');

*         bshunt.l(s,k)= sdata(s,'bsh0') ;
*         bshunt.lo(s,k)= sdata(s,'bsh0') ;
*         bshunt.up(s,k)= sdata(s,'bsh0') ;

         ratio.l(t,k)= tdata(t,'ratio') ;
         ratio.lo(t,k)= tdata(t,'ratio') ;
         ratio.up(t,k)= tdata(t,'ratio') ;

Equations

         ActPwrBalance_fn(i,k) Active power balance equation
         ReactPwrBalance_fn(i,k) Reactive power balance equation

*         ActPwrFlow_fn(f,i,k) Active power flow
*         ReactPwrFlow_fn(f,i,k) Active power flow

         GenLimit_P_L_fn(g,k) Limit on active DG generation (<)
         GenLimit_P_G_fn(g,k) Limit on active DG generation (>)
         GenLimit_Q_L_fn(g,k) Limit on reactive DG generation (<)
         GenLimit_Q_G_fn(g,k) Limit on reactive DG generation (>)
         FixPg_fn(g,k) Fix generator active power 
         FixQg_fn(g,k) Fix generator reactive power 

*         Current_fn(f,i,j,k) Longitudinal current limit 

*         VoltageLimit_G_fn(i,k) Voltage limit (>)
*         VoltageLimit_L_fn(i,k) Voltage limit (<)

*         GenVoltage_fn(g,i,k) imposed voltage at (substation) generators

*         Total_lines_connected_fn(k) 
*         Fix_disconnected_lines_fn(f,k)
*         Fix_line_status_fn(f,k) Fix the status of line breaker to intial values

*         Bstatuschange1_fn(f,k) expresses the change in branch breaker status   
*         Bstatuschange2_fn(f,k) expresses the change in branch breaker status 
*         Bstatus_lim_L(f,k) 
*         Bstatus_lim_U(f,k) 

*         Computevoltages(i,k) 

         Cost_fn Objective function

         ;

ActPwrBalance_fn(i,k)..
         sum(g $ (gendata(g,'Bus') eq idata(i,'BusN')), Pg(g,k)            ) -
         sum(l $ (ldata(l,'Bus')   eq idata(i,'BusN')), ldata(l,'PLoad')) -
         sum((f,j) $ ((bdata(f,'To')    eq idata(i,'BusN') and (bdata(f,'From') eq idata(j,'BusN')))), bdata(f,'BStatus0')*((e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k)-e(i,k)*e(j,k)-Vim(i,k)*Vim(j,k))*bdata(f,'g')+(e(i,k)*Vim(j,k)-Vim(i,k)*e(j,k))*bdata(f,'b'))) -
         sum((f,j) $ ((bdata(f,'From')  eq idata(i,'BusN') and (bdata(f,'To')   eq idata(j,'BusN')))), bdata(f,'BStatus0')*((e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k)-e(i,k)*e(j,k)-Vim(i,k)*Vim(j,k))*bdata(f,'g')+(e(i,k)*Vim(j,k)-Vim(i,k)*e(j,k))*bdata(f,'b'))) -
         sum((t,j) $ ((tdata(t,'To')    eq idata(i,'BusN') and (tdata(t,'From') eq idata(j,'BusN')))), (e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*tdata(t,'g')/(ratio(t,k)*ratio(t,k))-((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'g')+(-e(i,k)*Vim(j,k)+e(j,k)*Vim(i,k))*tdata(t,'b'))/ratio(t,k)) -
         sum((t,j) $ ((tdata(t,'From')  eq idata(i,'BusN') and (tdata(t,'To')   eq idata(j,'BusN')))), (e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*tdata(t,'g')-((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'g')+(-e(i,k)*Vim(j,k)+e(j,k)*Vim(i,k))*tdata(t,'b'))/ratio(t,k)) =e= 0 ;

*         sum((t,j) $ ((tdata(t,'To')    eq idata(i,'BusN') and (tdata(t,'From') eq idata(j,'BusN')))), (e(j,k)*e(j,k)+Vim(j,k)*Vim(j,k))*tdata(t,'g')/(ratio(t,k)*ratio(t,k))-((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'g')+(-e(i,k)*Vim(j,k)+e(j,k)*Vim(i,k))*tdata(t,'b'))/ratio(t,k)) -
*         sum((t,j) $ ((tdata(t,'From')  eq idata(i,'BusN') and (tdata(t,'To')   eq idata(j,'BusN')))), (e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*tdata(t,'g')-((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'g')+(-e(i,k)*Vim(j,k)+e(j,k)*Vim(i,k))*tdata(t,'b'))/ratio(t,k)) =e= 0 ;

ReactPwrBalance_fn(i,k)..
         sum(g $ (gendata(g,'Bus') eq idata(i,'BusN')), Qg(g,k)            ) -
         sum(l $ (ldata(l,'Bus')   eq idata(i,'BusN')), ldata(l,'QLoad')) +  
*         sum(s $ (sdata(s,'Bus')   eq idata(i,'BusN')), bshunt(s,k)*(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))) -
         sum((f,j) $ ((bdata(f,'To')    eq idata(i,'BusN') and (bdata(f,'From') eq idata(j,'BusN')))), bdata(f,'BStatus0')*(-(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*(bdata(f,'b')+bdata(f,'bsh'))+(e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*bdata(f,'b')+(e(i,k)*Vim(j,k)-Vim(i,k)*e(j,k))*bdata(f,'g'))) -
         sum((f,j) $ ((bdata(f,'From')  eq idata(i,'BusN') and (bdata(f,'To')   eq idata(j,'BusN')))), bdata(f,'BStatus0')*(-(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*(bdata(f,'b')+bdata(f,'bsh'))+(e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*bdata(f,'b')+(e(i,k)*Vim(j,k)-Vim(i,k)*e(j,k))*bdata(f,'g'))) -
         sum((t,j) $ ((tdata(t,'To')    eq idata(i,'BusN') and (tdata(t,'From') eq idata(j,'BusN')))), -(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*tdata(t,'b')/(ratio(t,k)*ratio(t,k))+((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'b')+(e(i,k)*Vim(j,k)-e(j,k)*Vim(i,k))*tdata(t,'g'))/ratio(t,k)) -
         sum((t,j) $ ((tdata(t,'From')  eq idata(i,'BusN') and (tdata(t,'To')   eq idata(j,'BusN')))), -(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*(tdata(t,'b')+tdata(t,'bsh'))+((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'b')+(e(i,k)*Vim(j,k)-e(j,k)*Vim(i,k))*tdata(t,'g'))/ratio(t,k)) =e= 0 ;

*         sum((t,j) $ ((tdata(t,'To')    eq idata(i,'BusN') and (tdata(t,'From') eq idata(j,'BusN')))), -(e(j,k)*e(j,k)+Vim(j,k)*Vim(j,k))*tdata(t,'b')/(ratio(t,k)*ratio(t,k))+((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'b')+(-e(i,k)*Vim(j,k)+e(j,k)*Vim(i,k))*tdata(t,'g'))/ratio(t,k)) -
*         sum((t,j) $ ((tdata(t,'From')  eq idata(i,'BusN') and (tdata(t,'To')   eq idata(j,'BusN')))), -(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k))*(tdata(t,'b')+tdata(t,'bsh'))+((e(i,k)*e(j,k)+Vim(i,k)*Vim(j,k))*tdata(t,'b')+(e(i,k)*Vim(j,k)-e(j,k)*Vim(i,k))*tdata(t,'g'))/ratio(t,k)) =e= 0 ;

GenLimit_P_L_fn(g,k) $ (gendata(g,'genP') eq 1) .. Pg(g,k) =L= gendata(g,'Pmax') ;
GenLimit_P_G_fn(g,k) $ (gendata(g,'genP') eq 1) .. Pg(g,k) =G= gendata(g,'Pmin') ;
GenLimit_Q_L_fn(g,k) $ (gendata(g,'genQ') eq 1) .. Qg(g,k) =L= gendata(g,'Qmax') ;
GenLimit_Q_G_fn(g,k) $ (gendata(g,'genQ') eq 1) .. Qg(g,k) =G= gendata(g,'Qmin') ;

FixPg_fn(g,k) $ (gendata(g,'genP') eq 0) .. Pg(g,k) =E= gendata(g,'P0') ;
FixQg_fn(g,k) $ (gendata(g,'genQ') eq 0) .. Qg(g,k) =E= gendata(g,'Q0') ;

*Current_fn(f,i,j,k) $ ( (idata(i,'BusN') EQ bdata(f,'From')) and (idata(j,'BusN') EQ bdata(f,'To'  )))..
*         (power(bdata(f,'g'),2) + power(bdata(f,'b'),2)) * (power(e(i,k)-e(j,k),2)+power(Vim(i,k)-Vim(j,k),2)) - (1-BStatus(f,p))*bigM =L= BStatus(f,p)*power(bdata(f,'Snom'),2);

*Current_fn(f,i,j,k) $((bdata(f,'Bstatus0') eq 1) and (idata(i,'BusN') EQ bdata(f,'From')) and (idata(j,'BusN') EQ bdata(f,'To'  )))..
*         (power(bdata(f,'g'),2) + power(bdata(f,'b'),2)) * (power(e(i,k)-e(j,k),2)+power(Vim(i,k)-Vim(j,k),2)) =L= power(bdata(f,'Snom'),2);

*VoltageLimit_G_fn(i,k).. e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k) =G= (1-VR)*idata(i,'Vmin')**2 ;
*VoltageLimit_L_fn(i,k).. e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k) =L= (1+VR)*idata(i,'Vmax')**2 ;

*GenVoltage_fn(g,i,k) $ ((gendata(g,'Bus') EQ  idata(i,'BusN'))) .. U(i,k) =E= power(gendata(g,'Vgen'),2) ;
*GenVoltage_fn(g,i,k) $ ((gendata(g,'Bus') EQ 0)) .. U(i,k) =E= power(gendata(g,'Vgen'),2) ;

* (longitudinal) current constraints 

*Total_lines_connected_fn(k).. sum(f,BStatus(f,k)) =E= NL ;

*Bstatuschange1_fn(f,k) .. BStatuschange(f,k) =G= BStatus(f,k) - bdata(f,'BStatus0');
*Bstatuschange2_fn(f,k) .. BStatuschange(f,k) =G= bdata(f,'BStatus0') - BStatus(f,k);

*Bstatus_lim_L(f,k) .. BStatus(f,k) =L= 1 ;
*Bstatus_lim_U(f,k) .. BStatus(f,k) =G= 0 ;

*Fix_disconnected_lines_fn(f,k) $ (bdata(f,'dBstatus') eq 1) .. BStatus(f,k) =e= bdata(f,'BStatus0') ;

*Fix_line_status_fn(f,k) .. BStatus(f,k) =e= bdata(f,'BStatus0');

*lim_switchings_fn(k) .. sum(f $ (bdata(f,'dBstatus') eq 1), BStatuschange(f,k)) =L= 2 ;

*Computevoltages(i,k) .. voltage(i,k) =E= sqrt(e(i,k)*e(i,k)+Vim(i,k)*Vim(i,k)) ;

*Cost_fn.. Cost =e= sum((g,k) $ (gendata(g,'genP') eq 1), 10000*power(Pg(g,k)-gendata(g,'P0'),2));

Cost_fn.. Cost =e= sum((g,k) $ (gendata(g,'genP') eq 1), gendata(g,'c0')+gendata(g,'c1')*100*Pg(g,k)+gendata(g,'c2')*10000*power(Pg(g,k),2)) ; 

model OPF_RC / all / ;

*restore_supply.optfile=1;

option profile = 3 ;
*option savepoint = 2;
*option minlp = sbb;
option nlp = conopt;
*option minlp = bonmin;
*option nodlim = 2000;
*OPF_RC.nodlim = 5000;
*OPF_RC.optcr = 0.;
solve OPF_RC using nlp minimizing cost ;
*solve OPF_RC using minlp minimizing cost ;
*solve OPF_RC using miqcp minimizing cost ;

*execute_unload "results.gdx" x.L x.M



